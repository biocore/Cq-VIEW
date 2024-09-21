import os
import pandas
import pathlib
from src.shared_symbols import *
from src.utils import sort_df_by_date, get_site_group_dict_by_name, \
    download_inputs_urls, extract_config_dict, get_input_and_output_dirs

PROTOCOL_NAME_KEY = "protocol_name"

STANDARD_FNAME_BASE = "standard_data"
# standardized output keys
EXTERNAL_SAMPLE_NAME_KEY = "external_sample_name"
EXTERNAL_SAMPLER_NAME_KEY = "external_sampler_name"
SAMPLE_COLLECTION_DATETIME_KEY = "sample_collection_datetime"
EXTERNAL_RESULT_NAME_KEY = "external_result_name"
SAMPLE_QUANT_DATETIME_KEY = "sample_quant_datetime"
EXTERNAL_PROTOCOL_NAME_KEY = "external_quant_protocol_name"
EXTERNAL_TARGET_NAME_KEY = "external_target_name"
TARGET_CQ_KEY = "target_cq"
TARGET_TO_CONTROL_CONC_RATIO_KEY = "target_to_control_conc_ratio"
TARGET_COPIES_PER_LITER_KEY = "target_copies_per_liter"
ACCEPT_KEY = "accept"

# campus dashboard output keys
CAMPUS_SAMPLE_KEY = EXTERNAL_SAMPLE_NAME_KEY
CAMPUS_SAMPLER_KEY = EXTERNAL_SAMPLER_NAME_KEY
CAMPUS_COLLECTION_DATETIME_KEY = SAMPLE_COLLECTION_DATETIME_KEY
CAMPUS_QPCR_RUN_NAME_KEY = EXTERNAL_RESULT_NAME_KEY
CAMPUS_QPCR_DATETIME_KEY = "sample_qpcr_datetime"
CAMPUS_QPCR_PROTOCOL_NAME_KEY = "external_qpcr_protocol_name"
CAMPUS_CQ_KEY = "sarscov2_cq"
CAMPUS_RATIO_KEY = "sarscov2_to_fecal_control_conc_ratio"
CAMPUS_CONC_KEY = "sarscov2_copies_per_liter"

# ddPCR-specific keys
DDPCR_WELL_KEY = "Well"
DDPCR_SAMPLE_DESC_1_KEY = "Sample description 1"
DDPCR_SAMPLE_DESC_2_KEY = "Sample description 2"
DDPCR_SAMPLE_DESC_3_KEY = "Sample description 3"
DDPCR_SAMPLE_DESC_4_KEY = "Sample description 4"
DDPCR_TARGET_KEY = "Target"
DDPCR_COPIES_PER_UL_KEY = "Conc(copies/µL)"
DDPCR_STATUS_KEY = "Status"
DDPCR_BAD_FLAG = "CHECK"

ORIGINAL_PREFIX = "original_"
REVISED_PREFIX = "revised_"
NORM_CONC_KEY = "norm_conc_key"
SITE_LOCATION_KEY = "site_location"

SCRIPT_TEMPLATE = \
    """
#!/usr/bin/env bash

# NB: before using, you must configure your github credentials by running
# gh auth login
#
# Also, I got permission denied fails on the gh pr create until I chmod'd
# read permissions for everyone (sudo chmod 644) on the
# config.yml and hosts.yml files in the <my_user>/.config/gh directory

# Uses variables:
# REPORT_DIR
# RUN_NAME
# REPO_DIR
# ANACONDADIR

# if any errors happened, the results shouldn't be pushed to the repository
if [ -s REPORT_DIR/generate_reports_error.log ]; then
  echo "Errors were logged while RUN_NAME processing; repo upload cancelled"
  exit 1
fi

# Add updated qpcr data files to local fork, then make PR for andersen lab repo
cd REPOS_DIR/SARS-CoV-2_WasteWater_San-Diego || exit
git checkout master
# get any updates from remote forked repo into local fork
git pull
# get any updates from remote andersen lab (original) repo into local original
git fetch andersen_wastewater
# merge any changes from the local original into the local fork
git merge --no-edit andersen_wastewater/master

cp REPORT_DIR/output/*_sewage_qPCR.csv \
    REPOS_DIR/SARS-CoV-2_WasteWater_San-Diego

# NB: no need to `git add` because we are only updating existing files
git commit -a -m "RUN_NAME"  # orig repo updates + local changes, to fork
git push

# Make a PR to the andersen lab (original) repo
source ANACONDADIR/activate wastewater
gh pr create --repo andersen-lab/SARS-CoV-2_WasteWater_San-Diego \
  --title "RUN_NAME" --body "Automated PR from Cq-VIEW"
source ANACONDADIR/deactivate

# TODO: copy the report dir into the deliveries dir?
"""


# write out to test campus dashboard file
def output_campus_df(curr_per_sample_all_w_viral_loads_df,
                     target_name, output_dir, timestamp_str, config):

    # NB: all dictionaries are ordered in Python 3.7+, so order of the columns
    # in the output file will be the same as the order of the keys in this dict
    cols_dict = {INTERNAL_SAMPLE_KEY: CAMPUS_SAMPLE_KEY,
                 LOCATION_KEY: CAMPUS_SAMPLER_KEY,
                 INTERNAL_COLLECTION_DATE_KEY: CAMPUS_COLLECTION_DATETIME_KEY,
                 RUN_DATE_KEY: CAMPUS_QPCR_DATETIME_KEY,
                 RESULT_ID_KEY: CAMPUS_QPCR_RUN_NAME_KEY,
                 CAMPUS_QPCR_PROTOCOL_NAME_KEY: CAMPUS_QPCR_PROTOCOL_NAME_KEY,
                 target_name: CAMPUS_CQ_KEY,
                 TARGET_TO_CONTROL_RATIO_KEY: CAMPUS_RATIO_KEY,
                 INTERNAL_NORM_CONC_KEY: CAMPUS_CONC_KEY,
                 ACCEPT_KEY: ACCEPT_KEY}

    base_out_df = _generate_base_qpcr_output_df(
        curr_per_sample_all_w_viral_loads_df, target_name, config)
    _trim_and_write_out_df(base_out_df, cols_dict, output_dir, timestamp_str,
                           "dashboard_data")


def output_standard_df_from_qpcr(
        curr_per_sample_all_w_viral_loads_df, assay_name, target_name,
        output_dir, timestamp_str, config):

    # NB: all dictionaries are ordered in Python 3.7+, so order of the columns
    # in the output file will be the same as the order of the keys in this dict
    cols_dict = {INTERNAL_SAMPLE_KEY: EXTERNAL_SAMPLE_NAME_KEY,
                 LOCATION_KEY: EXTERNAL_SAMPLER_NAME_KEY,
                 INTERNAL_COLLECTION_DATE_KEY: SAMPLE_COLLECTION_DATETIME_KEY,
                 RUN_DATE_KEY: SAMPLE_QUANT_DATETIME_KEY,
                 RESULT_ID_KEY: EXTERNAL_RESULT_NAME_KEY,
                 CAMPUS_QPCR_PROTOCOL_NAME_KEY: EXTERNAL_PROTOCOL_NAME_KEY,
                 EXTERNAL_TARGET_NAME_KEY: EXTERNAL_TARGET_NAME_KEY,
                 target_name: TARGET_CQ_KEY,
                 TARGET_TO_CONTROL_RATIO_KEY: TARGET_TO_CONTROL_CONC_RATIO_KEY,
                 INTERNAL_NORM_CONC_KEY: TARGET_COPIES_PER_LITER_KEY,
                 ACCEPT_KEY: ACCEPT_KEY}

    base_out_df = _generate_base_qpcr_output_df(
        curr_per_sample_all_w_viral_loads_df, target_name, config)
    base_out_df[EXTERNAL_TARGET_NAME_KEY] = assay_name

    _trim_and_write_out_df(base_out_df, cols_dict,
                           output_dir, timestamp_str, STANDARD_FNAME_BASE)


def output_standard_df_from_ddpcr(
        ddpcr_raw_out_df, control_name, protocol_name, protocol_date,
        output_dir, timestamp_str, out_name_base):

    # NB: all dictionaries are ordered in Python 3.7+, so order of the columns
    # in the output file will be the same as the order of the keys in this dict
    cols_dict = {EXTERNAL_SAMPLE_NAME_KEY: EXTERNAL_SAMPLE_NAME_KEY,
                 EXTERNAL_SAMPLER_NAME_KEY: EXTERNAL_SAMPLER_NAME_KEY,
                 SAMPLE_COLLECTION_DATETIME_KEY:
                     SAMPLE_COLLECTION_DATETIME_KEY,
                 SAMPLE_QUANT_DATETIME_KEY: SAMPLE_QUANT_DATETIME_KEY,
                 EXTERNAL_RESULT_NAME_KEY: EXTERNAL_RESULT_NAME_KEY,
                 EXTERNAL_PROTOCOL_NAME_KEY: EXTERNAL_PROTOCOL_NAME_KEY,
                 EXTERNAL_TARGET_NAME_KEY: EXTERNAL_TARGET_NAME_KEY,
                 TARGET_CQ_KEY: TARGET_CQ_KEY,
                 TARGET_TO_CONTROL_CONC_RATIO_KEY:
                     TARGET_TO_CONTROL_CONC_RATIO_KEY,
                 TARGET_COPIES_PER_LITER_KEY: TARGET_COPIES_PER_LITER_KEY,
                 ACCEPT_KEY: ACCEPT_KEY}

    standardized_df = _parse_raw_ddpcr_to_standardized_df(
        ddpcr_raw_out_df, control_name)
    standardized_df[EXTERNAL_PROTOCOL_NAME_KEY] = protocol_name
    standardized_df[SAMPLE_QUANT_DATETIME_KEY] = protocol_date
    standardized_df[TARGET_CQ_KEY] = ""

    out_base = f"{out_name_base}_ddPCR_{STANDARD_FNAME_BASE}"
    _trim_and_write_out_df(standardized_df, cols_dict,
                           output_dir, timestamp_str, out_base)


def update_location_viral_load_files(
        sample_base_df, site_group, source_dir, config_fp=None):

    config = extract_config_dict(config_fp)
    site_group_dict = get_site_group_dict_by_name(site_group, config)
    filename_pattern = site_group_dict[FILENAME_PATTERN_KEY]
    if SITE_LOCATIONS_KEY not in site_group_dict:
        raise ValueError(f"No site locations for site group "
                         f"'{site_group_dict[SITE_GROUP_KEY]}'")
    site_locations_by_prefix = site_group_dict[SITE_LOCATIONS_KEY]
    collection_date_key = config[OUTPUT_FILES_KEY][COLLECTION_DATE_KEY]

    input_dir, output_dir = get_input_and_output_dirs(source_dir)
    download_inputs_urls(input_dir)

    prev_report_dfs_by_report_name = _load_past_location_viral_load_dfs(
        input_dir, filename_pattern, site_locations_by_prefix)

    extended_fps = _extend_location_results_files(
        prev_report_dfs_by_report_name, output_dir, sample_base_df,
        filename_pattern, site_locations_by_prefix, collection_date_key,
        config)

    return extended_fps


def make_repo_upload_script(report_name, report_dir, repo_dir, conda_bin_dir):
    script_text = SCRIPT_TEMPLATE
    script_text = script_text.replace("RUN_NAME", report_name)
    script_text = script_text.replace("REPORT_DIR", report_dir)
    script_text = script_text.replace("REPOS_DIR", repo_dir)
    script_text = script_text.replace("ANACONDADIR", conda_bin_dir)

    _, output_dir = get_input_and_output_dirs(report_dir)
    output_fp = os.path.join(output_dir, f"{report_name}_repo_upload.sh")
    with open(output_fp, "w") as fh:
        fh.write(script_text)


def _generate_base_qpcr_output_df(curr_per_sample_all_w_viral_loads_df,
                                  target_name, config):

    working_df = curr_per_sample_all_w_viral_loads_df.copy()
    working_df[CAMPUS_QPCR_PROTOCOL_NAME_KEY] = config[PROTOCOL_NAME_KEY]
    working_df[ACCEPT_KEY] = ~working_df[IGNORE_KEY]

    # Don't let PEP change these to "if X" instead of "if X == False" since
    # the former doesn't work for pandas columns
    is_not_ignore_mask = working_df[IGNORE_KEY] == False  # noqa E712
    is_not_target_present_mask = \
        working_df[TARGET_PRESENT_KEY] == False  # noqa E712
    show_neg_one_mask = is_not_ignore_mask & is_not_target_present_mask
    working_df.loc[show_neg_one_mask, target_name] = \
        INTERNAL_MEASURED_BUT_NOT_FOUND_VALUE

    return working_df


def _trim_and_write_out_df(output_df, cols_dict,
                           output_dir, timestamp_str, out_base_str):
    trimmed_df = output_df.loc[:, cols_dict.keys()].copy()
    trimmed_df.rename(columns=cols_dict, inplace=True)
    trimmed_df.to_csv(
        f"{output_dir}/{timestamp_str}_{out_base_str}.csv", index=False)


def _load_past_location_viral_load_dfs(
        input_dir, filename_pattern, site_locations_by_prefix):
    prev_report_df_by_report_name = {}
    for curr_prefix, curr_location in site_locations_by_prefix.items():
        curr_report_fname = filename_pattern.format(curr_location)
        curr_prev_report_fp = f"{input_dir}/{curr_report_fname}"
        curr_prev_report_df = pandas.read_csv(curr_prev_report_fp)
        prev_report_df_by_report_name[curr_prefix] = curr_prev_report_df
    # next location

    return prev_report_df_by_report_name


def _extend_location_results_files(
        prev_report_df_by_report_name, output_dir, sample_base_df,
        filename_pattern, site_locations_by_prefix, collection_date_key,
        config):
    output_fps = []
    output_path = pathlib.Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    for curr_prefix, curr_location in site_locations_by_prefix.items():
        new_df = _prep_output_df(sample_base_df, config, [curr_prefix])

        curr_report_fname = filename_pattern.format(curr_location)
        curr_output_fp = f"{output_dir}/{curr_report_fname}"

        curr_prev_report_df = prev_report_df_by_report_name[curr_prefix]

        date_indexed_new_df = new_df.set_index(collection_date_key, drop=False)
        date_indexed_curr_prev_report_df = curr_prev_report_df.set_index(
            collection_date_key, drop=False)

        # use combine_first so that any reruns of a given sample date in
        # the new df will replace that sample date's row from the old dataframe
        output_df = date_indexed_new_df.combine_first(
            date_indexed_curr_prev_report_df)
        # remove any records that for any reason don't have a collection date
        # ... like the sample ENCUNKNOWN <smh ...>
        output_df = output_df[output_df[collection_date_key].notna()].copy()
        output_df = _polish_df_format(output_df, collection_date_key)

        output_df.to_csv(curr_output_fp, index=False, float_format='%.0f')
        output_fps.append(curr_output_fp)
    # next location

    return output_fps


def _prep_output_df(sample_base_df, config, site_prefix_list):
    location_mask = sample_base_df[SAMPLE_BASE_KEY].str.startswith(
        tuple(site_prefix_list))
    not_na_mask = sample_base_df[INTERNAL_NORM_CONC_KEY].notna()
    not_ignore_mask = sample_base_df[IGNORE_KEY] == False  # noqa E712
    combined_mask = location_mask & not_na_mask & not_ignore_mask

    cols_to_output = [INTERNAL_COLLECTION_DATE_KEY, INTERNAL_NORM_CONC_KEY]
    location_df = sample_base_df.loc[combined_mask, cols_to_output]

    date_format = '%m/%d/%y'
    location_df[INTERNAL_COLLECTION_DATE_KEY] = \
        location_df[INTERNAL_COLLECTION_DATE_KEY].dt.strftime(date_format)

    collection_date_key = config[OUTPUT_FILES_KEY][COLLECTION_DATE_KEY]
    norm_conc_key = config[OUTPUT_FILES_KEY][NORM_CONC_KEY]
    location_df.rename(columns={
        INTERNAL_COLLECTION_DATE_KEY: collection_date_key,
        INTERNAL_NORM_CONC_KEY: norm_conc_key}, inplace=True)

    return location_df


def _polish_df_format(a_df, collection_date_key):
    prepped_df = a_df.copy().round(0)
    prepped_df.reset_index(inplace=True, drop=True)
    prepped_df = sort_df_by_date(prepped_df, collection_date_key)
    prepped_df.fillna(0, inplace=True)
    return prepped_df


def write_excel(processing_run_name, dirnames, base_dfs, sample_dfs, log_fps,
                output_dir):
    excel_writer = pandas.ExcelWriter(
        f"{output_dir}/{processing_run_name}.xlsx",
        engine='xlsxwriter',
        datetime_format='mm/dd/yy',
        date_format='mm/dd/yy'
    )

    # Get the xlsxwriter workbook and worksheet objects.
    # NB: Pycharm says this attribute reference is unresolved, but it really
    # *is* there ...
    workbook = excel_writer.book

    # Add some cell formats.

    workbook.formats[0].set_font_size(10)
    workbook.formats[0].set_font_name('Courier New')
    whole_num_w_separator = workbook.add_format(
        {'num_format': '#,##', 'font_name': 'Courier New', 'font_size': 10})

    # Note: It isn't possible to format any cells that already have a format
    # like the index or headers or any cells that contain dates or datetimes.

    for curr_idx in range(0, len(dirnames)):
        if base_dfs[curr_idx] is not None:
            trunc_dirname = dirnames[curr_idx][:25]
            _add_to_excel(base_dfs[curr_idx], excel_writer,
                          f"sb_{trunc_dirname}", whole_num_w_separator)

    for curr_idx in range(0, len(dirnames)):
        if sample_dfs[curr_idx] is not None:
            trunc_dirname = dirnames[curr_idx][:25]
            _add_to_excel(sample_dfs[curr_idx], excel_writer,
                          f"s_{trunc_dirname}", whole_num_w_separator)

    for curr_idx in range(0, len(dirnames)):
        if log_fps[curr_idx] is not None:
            trunc_dirname = dirnames[curr_idx][:25]
            log_fp = log_fps[curr_idx]
            with open(log_fp) as fh:
                log_lines = fh.readlines()
            log_lines_df = pandas.DataFrame(log_lines)
            _add_to_excel(log_lines_df, excel_writer, f"cfg_{trunc_dirname}",
                          None, header=False)

    excel_writer.save()


def _add_to_excel(a_df, excel_writer, sheet_name, conc_format, header=True):
    date_columns = a_df.select_dtypes(include=['datetime64[ns, UTC]']).columns
    for date_column in date_columns:
        a_df[date_column] = a_df[date_column].dt.date

    a_df.to_excel(excel_writer, sheet_name, index=False, header=header)

    if conc_format is not None:
        worksheet = excel_writer.sheets[sheet_name]
        worksheet.freeze_panes(1, 3)  # Freeze the first row and 2nd column
        worksheet.set_column(2, 2, 20)  # make the name column wider
        for curr_col_idx in range(0, len(a_df.columns)):
            curr_col_name = a_df.columns[curr_col_idx]
            if CONC_SUFFIX in curr_col_name or \
                    curr_col_name == INTERNAL_NORM_CONC_KEY:
                worksheet.set_column(
                    curr_col_idx, curr_col_idx, 20, conc_format)


def _parse_raw_ddpcr_to_standardized_df(
        ddpcr_raw_out_df, control_name,
        sample_name_key=DDPCR_SAMPLE_DESC_1_KEY,
        sampler_name_key=DDPCR_SAMPLE_DESC_2_KEY,
        sample_collected_date_key=DDPCR_SAMPLE_DESC_3_KEY,
        result_name_key=DDPCR_SAMPLE_DESC_4_KEY):

    results_dicts_list = []

    # group the data by sample in well
    grouped = ddpcr_raw_out_df.groupby([sample_name_key, DDPCR_WELL_KEY])
    for curr_well_name, _ in grouped:
        curr_well_df = grouped.get_group(curr_well_name)

        well_grouped = curr_well_df.groupby([DDPCR_TARGET_KEY])
        group_sizes = well_grouped.size()
        if group_sizes.max() > 1:
            raise ValueError(
                f"Expected only one row per sample per target per well but "
                f"found multiple for '{curr_well_name}'")

        # get the sample info for the group
        curr_well_sample_name = curr_well_df[sample_name_key].values[0]
        curr_well_sampler_name = curr_well_df[sampler_name_key].values[0]
        curr_well_sample_collected_date = \
            curr_well_df[sample_collected_date_key].values[0]
        curr_well_result_name = curr_well_df[result_name_key].values[0]

        # get the concentration of the control
        curr_control_accept, curr_control_conc = \
            _get_accept_and_concentration(curr_well_df, control_name)

        for curr_target_name_in_well in \
                curr_well_df[DDPCR_TARGET_KEY].unique():
            if curr_target_name_in_well == control_name:
                continue

            # get the concentration of the target and calc the ratio
            curr_target_accept, curr_target_conc = \
                _get_accept_and_concentration(
                    curr_well_df, curr_target_name_in_well)
            curr_target_to_control_ratio = \
                curr_target_conc/curr_control_conc if curr_control_accept \
                else ""

            # make a results dict
            curr_results_dict = {
                DDPCR_WELL_KEY: curr_well_name,
                EXTERNAL_SAMPLE_NAME_KEY: curr_well_sample_name,
                EXTERNAL_SAMPLER_NAME_KEY: curr_well_sampler_name,
                SAMPLE_COLLECTION_DATETIME_KEY:
                    curr_well_sample_collected_date,
                EXTERNAL_RESULT_NAME_KEY: curr_well_result_name,
                EXTERNAL_TARGET_NAME_KEY: curr_target_name_in_well,
                TARGET_TO_CONTROL_CONC_RATIO_KEY: curr_target_to_control_ratio,
                # note 10^6 is conversion factor from copies/ul to copies/L
                TARGET_COPIES_PER_LITER_KEY: curr_target_conc * (10**6),
                ACCEPT_KEY: curr_target_accept
            }
            results_dicts_list.append(curr_results_dict)
        # next target in well
    # next well

    return pandas.DataFrame(results_dicts_list)


def _get_accept_and_concentration(well_df, target_name):
    # get the status of the target
    target_status = well_df.loc[well_df[DDPCR_TARGET_KEY] == target_name,
            DDPCR_STATUS_KEY].values[0]
    if target_status == DDPCR_BAD_FLAG:
        return False, pandas.NA

    # get the concentration of the control
    target_conc = well_df.loc[well_df[DDPCR_TARGET_KEY] == target_name,
            DDPCR_COPIES_PER_UL_KEY].values[0]
    try:
        target_conc = float(target_conc)
    except ValueError:
        target_conc = pandas.NA

    return True, target_conc