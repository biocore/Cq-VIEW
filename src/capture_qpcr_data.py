import copy
import glob
import logging
import numpy
import os
import pandas
import pathlib
import re
from datetime import datetime, timedelta
from dateutil import parser
from dateutil.relativedelta import relativedelta
from src.shared_symbols import *
from src.utils import extract_yaml_dict


def get_cq_conf_key_for_target(target_name, config):
    cq_conf_key = config[INPUT_FILES_KEY][CONF_KEY]
    target_cq_conf_key = f"{cq_conf_key}_{target_name}"
    return target_cq_conf_key


def get_unique_sorted_date_tuples(cum_per_sample_df, sort_asc=True):
    # get unique dates in the dataset
    dataset_ids_df = cum_per_sample_df.loc[:, [RUN_NAME_KEY, RUN_DATE_KEY]]
    unique_dataset_ids_df = dataset_ids_df.drop_duplicates().copy()
    unique_dataset_ids_df.sort_values(by=[RUN_DATE_KEY], inplace=True,
                                      ascending=sort_asc)
    sorted_unique_dataset_id_tuples = list(unique_dataset_ids_df.itertuples(
        index=False))
    return sorted_unique_dataset_id_tuples


def get_sorted_qpcr_results_fps(inputs_dir, config,
                                config_key=FILENAME_PATTERN_KEY):
    qpcr_fname_pattern = config[INPUT_FILES_KEY][config_key]

    qpcr_results_fps = glob.glob(f"{inputs_dir}/{qpcr_fname_pattern}")
    qpcr_results_fps.sort(key=os.path.getmtime, reverse=False)
    return qpcr_results_fps


def get_cumulative_data(inputs_dir, output_dir, report_name, config):
    # load the data
    per_sample_df_by_date = capture_qpcr_data(inputs_dir, config)
    if not per_sample_df_by_date:
        logger = logging.getLogger()
        logger.warning(f"No qpcr data found in {inputs_dir}")
        return None

    # combine all the data into one cumulative dataframe
    cumulative_raw_per_sample_df = combine_all_datesets(per_sample_df_by_date)

    # output the cumulative raw data to a file
    os.makedirs(output_dir, exist_ok=True)
    cumulative_raw_per_sample_df.to_csv(
        f"{output_dir}/{report_name}_sample_data_complete.csv", index=False)
    return cumulative_raw_per_sample_df


def capture_qpcr_data(inputs_dir, config):
    per_sample_df_by_date = {}
    exclude_dicts_list = _load_excludes(inputs_dir, config)
    qpcr_results_fps = get_sorted_qpcr_results_fps(inputs_dir, config)
    for curr_qpcr_results_fp in qpcr_results_fps:
        logger = logging.getLogger()
        logger.info(f"parsing {curr_qpcr_results_fp}")

        curr_qpcr_run_name = pathlib.PurePath(curr_qpcr_results_fp).stem

        use_xlsx = curr_qpcr_results_fp.endswith("xlsx")
        curr_qpcr_run_date = _get_run_date_from_file_contents(
            curr_qpcr_results_fp, config, use_xlsx)
        curr_qpcr_df = _read_in_qpcr_output(
            curr_qpcr_results_fp, config, use_xlsx)

        curr_per_sample_df = _get_per_sample_df(
            curr_qpcr_df, curr_qpcr_run_name, curr_qpcr_run_date, config)
        curr_per_sample_df = _mark_excluded_records(
            curr_per_sample_df, exclude_dicts_list)
        per_sample_df_by_date[curr_qpcr_run_date] = curr_per_sample_df
    # next input qpcr file

    return per_sample_df_by_date


def _load_excludes(inputs_dir, config):
    exclude_dicts_list = []
    curr_excludes_fps = get_sorted_qpcr_results_fps(
        inputs_dir, config, config_key=EXCLUDE_FNAME_PATTERN_KEY)
    for curr_exclude_fp in curr_excludes_fps:
        curr_exclude_dict = extract_yaml_dict(curr_exclude_fp)
        exclude_dicts_list.append(curr_exclude_dict)

    return exclude_dicts_list


def _mark_excluded_records(per_sample_df, exclude_dicts_list):
    for curr_exclude_dict in exclude_dicts_list:
        curr_affected_run = curr_exclude_dict[RUN_NAME_KEY]
        curr_run_mask = per_sample_df[RUN_NAME_KEY] == curr_affected_run

        curr_excludes_list = curr_exclude_dict[EXCLUDE_KEY]
        for curr_exclude_details_dict in curr_excludes_list:
            curr_note = curr_exclude_details_dict[NOTE_KEY]
            for curr_key, curr_val in curr_exclude_details_dict.items():
                if curr_key != NOTE_KEY:
                    curr_exclude_mask = per_sample_df[curr_key] == curr_val
                    curr_combined_mask = curr_run_mask & curr_exclude_mask
                    per_sample_df.loc[curr_combined_mask, EXCLUDE_NOTE_KEY] = \
                        curr_note
                # endif key is not the note key
            # next key in current exclude dict
        # next exclude details dict in list
    # next exclude dict in list

    return per_sample_df


def _read_in_qpcr_output(qpcr_output_fp, config, use_xlsx=True):
    input_dict = config[INPUT_FILES_KEY]
    sample_key = input_dict[SAMPLE_KEY]
    target_key = input_dict[TARGET_KEY]

    if use_xlsx:
        # read in just the tabular data from the xlsx, skipping settings info
        qpcr_output_df = pandas.read_excel(
            io=qpcr_output_fp, sheet_name=input_dict[SHEET_NAME_KEY],
            skiprows=input_dict[ROWS_TO_SKIP_KEY])
    else:
        qpcr_output_df = pandas.read_csv(
            qpcr_output_fp, comment="#")

    qpcr_output_df.dropna(subset=[sample_key, target_key], inplace=True)

    return qpcr_output_df


def _get_run_date_from_file_contents(qpcr_output_fp, config, use_xlsx):
    date_created_key = config[INPUT_FILES_KEY][DATE_CREATED_KEY_NAME]
    rows_to_read = config[INPUT_FILES_KEY][ROWS_TO_SKIP_KEY] - 1
    metadata_separator = config[INPUT_FILES_KEY][METADATA_SEPARATOR_KEY]
    setting_name_key = "setting_name"
    setting_val_key = "setting_val"

    # read in the file from the beginning and dig the run date out
    # of the settings info at the top
    if use_xlsx:
        qpcr_output_df = pandas.read_excel(
            sheet_name=config[INPUT_FILES_KEY][SHEET_NAME_KEY],
            io=qpcr_output_fp, header=None, nrows=rows_to_read)
    else:
        qpcr_output_df = pandas.read_csv(
            qpcr_output_fp, header=None, sep=metadata_separator,
            nrows=rows_to_read)

    qpcr_output_df = qpcr_output_df.rename(
        columns={
            qpcr_output_df.columns[0]: setting_name_key,
            qpcr_output_df.columns[1]: setting_val_key}
    )

    # endswith instead of == because when reading .csv, the metadata lines all
    # start with "# " ...
    desired_row_mask = \
        qpcr_output_df[setting_name_key].str.endswith(date_created_key)
    run_date_series = qpcr_output_df.loc[desired_row_mask, setting_val_key]
    run_date_str = run_date_series.iloc[0]
    run_date = parser.parse(run_date_str, fuzzy=True)
    run_date = run_date.replace(tzinfo=None)

    return run_date


def _get_per_sample_df(qpcr_output_df, qpcr_run_name, qpcr_run_date, config):
    null_value = config[INPUT_FILES_KEY][NULL_VALUE_KEY]
    well_key = config[INPUT_FILES_KEY][WELL_KEY]
    sample_key = config[INPUT_FILES_KEY][SAMPLE_KEY]
    target_key = config[INPUT_FILES_KEY][TARGET_KEY]
    raw_cqs_key = config[INPUT_FILES_KEY][RAW_CQ_KEY]
    conf_key = config[INPUT_FILES_KEY][CONF_KEY]

    cleaned_df = qpcr_output_df.copy()
    cleaned_df.reset_index(drop=False, inplace=True)
    cleaned_df[INTERNAL_SAMPLE_KEY] = \
        cleaned_df[sample_key]
    cleaned_df[INTERNAL_WELL_KEY] = \
        cleaned_df[well_key]
    cleaned_df[INTERNAL_SAMPLE_AND_WELL_KEY] = \
        cleaned_df[INTERNAL_SAMPLE_KEY] + "." + \
        cleaned_df[INTERNAL_WELL_KEY] + "." + qpcr_run_name
    # TODO: take out this hack for past input typos or leave it in?
    cleaned_df[target_key] = \
        cleaned_df[target_key].str.replace('PPMoV', 'PMMoV')

    try:
        per_sample_df = cleaned_df.pivot(
            index=[INTERNAL_WELL_KEY, INTERNAL_SAMPLE_AND_WELL_KEY,
                   INTERNAL_SAMPLE_KEY],
            columns=target_key,
            values=[raw_cqs_key, conf_key])
    except ValueError:
        logger = logging.getLogger()
        logger.warning("Non-unique well/sample/target combinations found")
        logger.warning(cleaned_df[[well_key, INTERNAL_SAMPLE_AND_WELL_KEY,
                                   target_key]].value_counts())
        raise

    # the col names are tuples made up of raw_cqs_key, target_key or
    # conf_key, target_key; here rename them to be easier to work with
    new_col_names = ["_".join(x).replace(f"{raw_cqs_key}_", "")
                     for x in per_sample_df.columns]
    per_sample_df.columns = new_col_names
    per_sample_df.replace(
        null_value, INTERNAL_MEASURED_BUT_NOT_FOUND_VALUE, inplace=True)
    per_sample_df = per_sample_df.astype(float)
    per_sample_df.reset_index(inplace=True)  # also adds sample name as column

    per_sample_df = _add_group_base_date_and_location(per_sample_df, config)
    per_sample_df[RUN_NAME_KEY] = qpcr_run_name
    per_sample_df[RUN_DATE_KEY] = qpcr_run_date
    per_sample_df[EXCLUDE_NOTE_KEY] = None

    return per_sample_df


def _add_group_base_date_and_location(sample_df, config):
    sample_info = sample_df.apply(
        lambda row: _extract_group_base_date_and_location(
            row[INTERNAL_SAMPLE_KEY], config),
        axis=1)

    # split column of lists into new columns then join back to orig df
    split = pandas.DataFrame(sample_info.to_list(),
                             columns=[SITE_GROUP_KEY, SAMPLE_BASE_KEY,
                                      LOCATION_KEY,
                                      INTERNAL_COLLECTION_DATE_KEY])
    sample_df = pandas.concat([sample_df, split], axis=1)

    return sample_df


def _extract_group_base_date_and_location(sample_id, config):
    found_site_group = found_site_prefix = None

    first_letter_and_after_regex = r".*?([A-z].*)$"
    re_match = re.match(first_letter_and_after_regex, sample_id)
    first_cap_and_after_str = re_match.group(1)

    site_groups = _get_site_groups(config)
    for curr_site_group_dict in site_groups:
        curr_site_group = curr_site_group_dict[SITE_GROUP_KEY]
        curr_site_prefixes = curr_site_group_dict[SITE_PREFIXES_KEY]
        for curr_site_prefix in curr_site_prefixes:
            if first_cap_and_after_str.startswith(curr_site_prefix):
                found_site_group = curr_site_group
                found_site_prefix = curr_site_prefix
                break
            # endif
        # next site prefix

        if found_site_group is not None:
            break
    # next group

    # TODO: figure out if these strings can be centralized somehow?
    if found_site_group == "sandiego":
        parse_method = _extract_sandiego_sample_base_location_and_date
    elif found_site_group in ["ucsd", "usd"]:
        parse_method = _extract_ucsd_sample_base_location_and_date
    elif found_site_group == "ucsf":
        parse_method = _extract_ucsf_sample_base_location_and_date
    elif found_site_group in [CONTROLS_GROUP, LADDERS_GROUP]:
        parse_method = _fake_ladder_sample_base_location_and_date
    elif found_site_group is None:
        if sample_id in ["unlabeled"]:
            parse_method = _fake_ignore_sample_base_location_and_date
        else:
            raise ValueError(
                f"No site group found for sample id '{sample_id}'")
    else:
        raise ValueError(f"Unrecognized site group '{found_site_group}'")

    sample_base, location, collection_date = parse_method(
        sample_id, found_site_prefix)

    return [found_site_group, sample_base, location, collection_date]


def _get_site_groups(config):
    ladder_prefixes = []
    for curr_assay_dict in config[INPUT_FILES_KEY][ASSAYS_KEY]:
        curr_targets_dict = curr_assay_dict[TARGETS_KEY]
        for curr_target_dict in curr_targets_dict:
            if TARGET_LADDER_PREFIX_KEY in curr_target_dict:
                ladder_prefixes.extend(
                    curr_target_dict[TARGET_LADDER_PREFIX_KEY])

    site_groups = copy.deepcopy(config[SITE_GROUPS_KEY])
    site_groups.extend([{
            SITE_GROUP_KEY: CONTROLS_GROUP,
            SITE_PREFIXES_KEY: config[INPUT_FILES_KEY][CONTROL_PREFIXES_KEY]
        },
        {
            SITE_GROUP_KEY: LADDERS_GROUP,
            SITE_PREFIXES_KEY: ladder_prefixes
        }]
    )

    return site_groups


def _extract_sandiego_sample_base_location_and_date(
        sample_id, site_prefix):
    # example sample base:
    # PLNOV29
    sample_id_pieces = sample_id.split(".")
    sample_base = sample_id_pieces[3]
    location = site_prefix
    date_base = sample_base.replace(site_prefix, "")
    collection_date = _extract_date_from_date_base(date_base)
    return sample_base, location, collection_date


# NB: site_prefix parameter is not used but is included to match signature
# of other extract functions used in same strategy pattern
def _extract_ucsd_sample_base_location_and_date(sample_id, site_prefix):
    # example sample base:
    # AS095
    # NB: this approach won't work for some very old ids (e.g., 12.8.AS015_V2)
    # from back before they started adding the year to the dates and when they
    # used to separate additional info from the sampler name with an underscore
    # instead of a dot.  However, I don't think they did a viral-load-capable
    # assay on any of these ids, so they shouldn't be coming through this
    # system anyway.
    sample_id_pieces = sample_id.split(".")
    sample_base = ".".join(sample_id_pieces[:4])  # nb upper bound is exclusive
    location = sample_id_pieces[3]

    date_base = "/".join(sample_id_pieces[:3])  # nb upper bound is exclusive
    prep_date = _extract_date_from_date_base(date_base)

    # Per RK 20230420: collection date should be day on which the bottle was
    # collected (and brought to the lab and prepped) not the day on which
    # the bottle was itself collecting liquid (i.e., the previous day, except
    # on weekends when it is really 3 days--fri, sat, sun).
    collection_date = prep_date
    # collection_date = prep_date - timedelta(days=1)
    return sample_base, location, collection_date


def _extract_ucsf_sample_base_location_and_date(sample_id, site_prefix):
    # example sample base:
    # 'OSP_12_19_22'
    sample_id_pieces = sample_id.split(".")
    # TODO: decide whether to put back. Older data have no replicates, so
    #  sample id is just, e.g., `OSP_10_31_22`, not `OSP_10_31_22.R1`
    # if len(sample_id_pieces) != 2:
    #     raise ValueError(f"ucsf sample_id '{sample_id}' splits into "
    #                      f"{len(sample_id_pieces)} pieces")
    sample_base = sample_id_pieces[0]
    location = site_prefix

    base_pieces = sample_base.split("_")
    if len(base_pieces) != 4:
        raise ValueError(f"ucsf sample base '{sample_base}' splits into "
                         f"{len(base_pieces)} pieces")

    date_base = "/".join(base_pieces[1:4])  # nb upper bound is exclusive
    collection_date = _extract_date_from_date_base(date_base)
    return sample_base, location, collection_date


def _fake_ladder_sample_base_location_and_date(sample_id, site_prefix):
    return sample_id, site_prefix, None


# NB: site_prefix parameter is not used but is included to match signature
# of other extract functions used in same strategy pattern
def _fake_ignore_sample_base_location_and_date(sample_id, site_prefix):
    return sample_id, None, None


def _extract_date_from_date_base(date_base):
    month_prefixes = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL",
                      "AUG", "SEP", "OCT", "NOV", "DEC"]
    try:
        if date_base[0].isalpha() and \
                not date_base.startswith(tuple(month_prefixes)):
            raise ValueError(f"Unknown date prefix in '{date_base}'")

        collection_date = parser.parse(date_base, fuzzy=True)
    except Exception: # noqa E722
        collection_date = None  # parser failures make this NaT
        logger = logging.getLogger()
        logger.warning(f"Cannot parse date base '{date_base}'")

    # Heuristic: deal with the fact that sample names don't have years in
    # them >:-( by assuming date in the future are actually from last yr.
    # See https://stackoverflow.com/a/51275729 for why not to use
    # timedelta (tldr: a yr isn't a fixed duration--some yrs are 366 days)
    if collection_date:
        current_datetime = datetime.now()
        if collection_date > current_datetime:
            collection_date = collection_date - relativedelta(years=1)

    return collection_date


def combine_all_datesets(per_sample_df_by_date):
    cumulative_samples_df = None

    for curr_date, curr_per_sample_df in per_sample_df_by_date.items():
        if cumulative_samples_df is None:
            cumulative_samples_df = curr_per_sample_df
        else:
            temp_concat_df = pandas.concat(
                [cumulative_samples_df, curr_per_sample_df],
                ignore_index=True)

            # TODO: decide whether to put back dedup that replaces earlier
            #  runs of a given sample with later runs
            # deduped_cum_df = temp_concat_df.sort_values(
            #     RUN_DATE_KEY, ascending=False).drop_duplicates(
            #     INTERNAL_SAMPLE_KEY).sort_index()
            # deduped_cum_df.reset_index(inplace=True, drop=True)
            # cum_sample_df = deduped_cum_df

            cumulative_samples_df = temp_concat_df.reset_index(drop=True)
        # end if
    # next per sample df

    return cumulative_samples_df


def extract_records_measured_for_target(input_df, target_name):
    # any records for which this target wasn't measured will be na;
    # exclude those!
    target_was_measured_mask = \
        input_df[target_name].notna()
    target_ladder_df = \
        input_df.loc[target_was_measured_mask, :].copy()
    return target_ladder_df


def convert_measured_but_not_found_to_nan(input_df):
    # now that the "not-measured" na's are gone,
    # replace any "measured but not found" values with na
    # so they are left out of calculations
    output_df = input_df.replace(
        INTERNAL_MEASURED_BUT_NOT_FOUND_VALUE, numpy.NaN)
    return output_df
