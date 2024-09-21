import logging
import os
import pandas
from datetime import datetime
from sys import argv
from src.shared_symbols import *
from src.utils import generate_timestamp_str, get_input_and_output_dirs, \
    extract_config_dict, set_up_logger, get_site_group_dict_by_name, \
    get_parent_dir
from src.relocate_local_files import centralize_run_files_from_dirs
from src.capture_and_plot_ladders import capture_and_save_ladders
from src.capture_qpcr_data import get_cumulative_data
from src.calculate_viral_loads import calc_viral_loads_for_site_group, \
    temp_ucsf_processing
from src.output_formatted_data import write_excel, \
    update_location_viral_load_files, make_repo_upload_script, \
    output_standard_df_from_ddpcr
from src.generate_plots import generate_plots

logging.basicConfig(level=logging.INFO)
ERRORLOG = logging.getLogger("ErrorLog")
ERRORLOG.setLevel(logging.ERROR)


def generate_ddpcr_report():
    """
    # example usage from command line:
    generate_ddpcr_report \
       "/Users/abirmingham/Work/Projects/covid_wastewater/ddpcr/reports/test" \
       "ddPCR" \
       "2022-09-13" \
       "PMMoV" \
       "" \
       "/Users/abirmingham/Desktop/220913_samples_data_A04_mod.csv"
    """

    # read in command-line arguments
    parent_output_dir = argv[1]
    protocol_name = argv[2]
    protocol_date = argv[3]
    control_name = argv[4]
    report_name_base = argv[5]
    source_fp = argv[6:]

    sources = {}
    for curr_source_fp in source_fp:
        # read curr_source_fp into a dataframe
        curr_source_df = pandas.read_csv(curr_source_fp)
        curr_source_fname = os.path.basename(curr_source_fp)
        curr_source_name = os.path.splitext(curr_source_fname)[0]
        sources[curr_source_name] = curr_source_df
    # next source file path

    ddpcr_raw_df = pandas.concat(sources.values(), ignore_index=True)

    if report_name_base == "":
        report_name_base = "_".join(sources.keys())
    timestamp_str = generate_timestamp_str()

    output_standard_df_from_ddpcr(ddpcr_raw_df, control_name, protocol_name,
                                  protocol_date, parent_output_dir,
                                  timestamp_str, report_name_base)


def generate_ucsf_report():
    """
    # example usage from command line:
    generate_ucsf_report \
       "/Applications/miniconda3/bin" \
       "/Users/abirmingham/Work/Repositories" \
       "/Users/abirmingham/Work/Projects/covid_wastewater/qpcr/data/production" \
       "/Users/abirmingham/Work/Projects/covid_wastewater/qpcr/reports" \
       "" \
       "230607_RTL_PL_ENC_SB"
    """

    # read in command-line arguments
    conda_bin_dir = argv[1]
    repo_dir = argv[2]
    parent_runs_dir = argv[3]
    parent_output_dir = argv[4]
    report_name_base = argv[5]
    earliest_date_str = None
    source_dir_names = argv[6:]

    parent_dir = get_parent_dir()
    config_fp = os.path.join(parent_dir, "config_ucsf.yml")

    _generate_report(conda_bin_dir, repo_dir, parent_runs_dir, parent_output_dir,
                     report_name_base, earliest_date_str, source_dir_names,
                     config_fp=config_fp)


def generate_report():
    """
    # example usage from command line:
    generate_report \
       "/Applications/miniconda3/bin" \
       "/Users/abirmingham/Work/Repositories" \
       "/Users/abirmingham/Work/Projects/covid_wastewater/qpcr/data/production" \
       "/Users/abirmingham/Work/Projects/covid_wastewater/qpcr/reports" \
       "" \
       "05/25/2023" \
       "230607_RTL_PL_ENC_SB"
    """

    # argv = ["",
    #         "/Applications/miniconda3/bin",
    #         "/Users/abirmingham/Work/Repositories",
    #         "/Users/abirmingham/Work/Projects/covid_wastewater/qpcr/data/production",
    #         "/Users/abirmingham/Work/Projects/covid_wastewater/qpcr/reports",
    #         "",
    #         "05/01/2023",
    #         "230602_UCSF"]

    # read in command-line arguments
    conda_bin_dir = argv[1]
    repo_dir = argv[2]
    parent_runs_dir = argv[3]
    parent_output_dir = argv[4]
    report_name_base = argv[5]
    earliest_date_str = argv[6]
    source_dir_names = argv[7:]

    _generate_report(conda_bin_dir, repo_dir, parent_runs_dir, parent_output_dir,
                     report_name_base, earliest_date_str, source_dir_names)


def _generate_report(conda_bin_dir, repo_dir, parent_runs_dir,
                     parent_output_dir, report_name_base, earliest_date_str,
                     source_dir_names, config_fp=None):
    report_name_base = "_".join(source_dir_names) if report_name_base == "" \
        else report_name_base
    timestamp_str = generate_timestamp_str()
    report_dir = os.path.join(
        parent_output_dir, f"{timestamp_str}_{report_name_base}")
    _, present_site_group_names = generate_source_dir_names_report(
        parent_runs_dir, report_dir, f"{report_name_base}_only",
        None, source_dir_names, repo_dir, conda_bin_dir, timestamp_str,
        config_fp=config_fp)

    config = extract_config_dict(config_fp)
    is_context_list = [_is_context_site_group(x, config) for x
                       in present_site_group_names]
    if any(is_context_list):
        generate_timeframe_report(
            parent_runs_dir, report_dir, f"{report_name_base}_cumulative",
            earliest_date_str, repo_dir, conda_bin_dir, timestamp_str)


def generate_timeframe_report(
        parent_runs_dir, parent_output_dir, report_name_base,
        earliest_date_str, repo_dir, conda_bin_dir, timestamp_str=None):

    earliest_date = _get_date_from_str(earliest_date_str)
    relevant_fps = _get_relevant_dir_paths(parent_runs_dir, earliest_date)

    timestamp_str, present_site_group_names = \
        _generate_report_from_source_dirs(
            relevant_fps, parent_output_dir, report_name_base,
            earliest_date, repo_dir, conda_bin_dir,
            timestamp_str=timestamp_str)
    return timestamp_str, present_site_group_names


def generate_source_dir_names_report(
        parent_runs_dir, parent_output_dir, report_name_base,
        earliest_date_str, source_runs_list, repo_dir, conda_bin_dir,
        timestamp_str=None, config_fp=None):
    earliest_date = _get_date_from_str(earliest_date_str)
    run_dirs = [os.path.join(parent_runs_dir, x) for x in source_runs_list]
    timestamp_str, present_site_group_names = \
        _generate_report_from_source_dirs(
            run_dirs, parent_output_dir, report_name_base,
            earliest_date, repo_dir, conda_bin_dir,
            timestamp_str=timestamp_str, config_fp=config_fp)
    return timestamp_str, present_site_group_names


def _get_date_from_str(date_str):
    if date_str is None or date_str.strip() == "":
        return None
    a_date = datetime.strptime(date_str, '%m/%d/%Y')
    return a_date


def _get_relevant_dir_paths(inputs_dir, earliest_date=None):
    # get all the folders that are immediate children of inputs_dir
    dir_paths = [os.path.join(inputs_dir, x) for x in os.listdir(inputs_dir)
                 if os.path.isdir(os.path.join(inputs_dir, x))]
    dir_paths.sort(key=os.path.getmtime, reverse=False)

    # include only paths modified on or after earliest_date
    if earliest_date is not None:
        temp_dir_paths = []
        for curr_dir_path in dir_paths:
            curr_date = datetime.fromtimestamp(os.path.getmtime(curr_dir_path))
            if curr_date >= earliest_date:
                temp_dir_paths.append(curr_dir_path)
        dir_paths = temp_dir_paths

    return dir_paths


def _generate_report_from_source_dirs(
        source_dirs_list, parent_output_dir, report_name_base,
        earliest_date, repo_dir, conda_bin_dir, timestamp_str=None,
        config_fp=None):
    try:
        timestamp_str = generate_timestamp_str() if timestamp_str is None \
            else timestamp_str
        report_name = f"{timestamp_str}_{report_name_base}"
        if earliest_date:
            earliest_date_str = earliest_date.strftime('%Y-%m-%d')
            report_name = f"{report_name}_since_{earliest_date_str}"

        report_dir = os.path.join(parent_output_dir, report_name)
        os.makedirs(report_dir)
        _add_error_log_filehandler(report_dir)

        centralize_run_files_from_dirs([None, source_dirs_list, report_dir])

        # TODO: Would be preferable to load the data once, then use for both
        #  capture_and_save_ladders and _report_viral_loads; current code
        #  loads individually in each

        # Heuristic: if there are more than about 3 ladder sets represented,
        # then adding the mixed model lines (along with the ladder lines and
        # the standard curve line) makes the graph too busy to read
        do_include_mixed_model = len(source_dirs_list) < 4
        capture_and_save_ladders(
            [report_dir], do_include_mixed_model, config_fp, report_name)

        # Limit the report to only site groups that require context if
        # an earliest date has been specified
        do_context_site_groups_only = earliest_date is not None
        present_site_group_names = _report_viral_loads(
            [report_dir, do_context_site_groups_only, report_name,
             repo_dir, conda_bin_dir, config_fp])

        return timestamp_str, present_site_group_names
    except Exception:  # noqa E722
        errorlog = logging.getLogger("ErrorLog")
        errorlog.exception("_generate_report_from_source_dirs failed")
        raise


def _add_error_log_filehandler(report_dir):
    handler = logging.FileHandler(
        os.path.join(report_dir, 'generate_reports_error.log'))
    ERRORLOG.addHandler(handler)


def _report_viral_loads(argslist):
    curr_source_fhandler = None
    logger = None
    try:
        # read in command-line arguments
        source_dir = argslist[0]
        do_context_site_groups_only = argslist[1]
        processing_run_name = argslist[2]
        repo_dir = argslist[3]
        conda_bin_dir = argslist[4]
        config_fp = argslist[5]

        # read in config file
        do_zoomed_plot = True
        config = extract_config_dict(config_fp)
        dirname = os.path.basename(source_dir)
        inputs_dir, output_dir = get_input_and_output_dirs(source_dir)

        log_fp = os.path.join(
            output_dir, f"{processing_run_name}_log_report_viral_loads.log")
        logger, curr_source_fhandler = set_up_logger(log_fp, config)
        logger.info(source_dir)

        cumulative_raw_per_sample_df = get_cumulative_data(
            inputs_dir, output_dir, processing_run_name, config)

        present_site_group_names = pandas.unique(
            cumulative_raw_per_sample_df[SITE_GROUP_KEY])
        for curr_site_group_name in present_site_group_names:
            if curr_site_group_name in [CONTROLS_GROUP, LADDERS_GROUP]:
                continue

            if curr_site_group_name == "ucsf" and config_fp is None:
                raise ValueError("Should not process ucsf data with default config!")
            if curr_site_group_name != "ucsf" and config_fp is not None:
                raise ValueError("Should not process non-ucsf data with non-default config!")

            if do_context_site_groups_only:
                is_context_site_group = \
                    _is_context_site_group(curr_site_group_name, config)
                if not is_context_site_group:
                    continue

            site_dir = os.path.join(output_dir, curr_site_group_name)
            _, site_output_dir = get_input_and_output_dirs(site_dir)
            sample_df, sample_base_df = _calc_for_site_group(
                curr_site_group_name, processing_run_name,
                cumulative_raw_per_sample_df, site_output_dir, logger,
                config_fp)

            if sample_df is not None:
                run_name_w_site = f"{processing_run_name}_" \
                                  f"{curr_site_group_name}"

                if curr_site_group_name == "ucsf":
                    temp_ucsf_processing(
                        sample_df, site_output_dir, run_name_w_site)
                else:

                    write_excel(
                        run_name_w_site, [dirname], [sample_base_df], [sample_df],
                        [log_fp], site_output_dir)

                    # TODO: make this more elegant. Only the sandiego site group
                    #  gets new results appended to running location files and
                    #  uploaded to github repos, and has per-location graphs made,
                    #  but maybe there's a better way to indicate this than a
                    #  hard-coded check ...
                    if curr_site_group_name == "sandiego":
                        update_location_viral_load_files(
                            sample_base_df, curr_site_group_name,
                            site_dir, config_fp)
                        make_repo_upload_script(
                            run_name_w_site, site_dir,
                            repo_dir, conda_bin_dir)

                        generate_plots(run_name_w_site, [dirname],
                                       [sample_df],
                                       [sample_base_df],
                                       site_output_dir, config,
                                       None, do_zoomed_plot)
                    # endif sandiego
                # endif curr_site_group_name == "ucsf" or not
            # endif calcs made
        # next site_group

        return present_site_group_names
    except Exception:
        logging.exception("_report_viral_loads failed")
        raise
    finally:
        if curr_source_fhandler:
            curr_source_fhandler.close()
            if logger:
                logger.removeHandler(curr_source_fhandler)


def _calc_for_site_group(curr_site_group_name, processing_run_name,
                         cumulative_raw_per_sample_df, site_output_dir,
                         logger, config_fp):
    config = extract_config_dict(config_fp)
    curr_site_group_dict = get_site_group_dict_by_name(
        curr_site_group_name, config)
    curr_report_name = \
        f"{processing_run_name}_{curr_site_group_name}"

    site_limited_per_sample_df = _extract_site_group_samples_only(
        cumulative_raw_per_sample_df, curr_site_group_name)
    if len(site_limited_per_sample_df) == 0:
        logger.info(f"No samples for site group '{curr_site_group_name}' "
                    f"in cumulative qPCR dataset")
        return None, None, None

    os.makedirs(site_output_dir, exist_ok=True)

    curr_per_sample_viral_loads_df, curr_cum_per_sample_base_df = \
        calc_viral_loads_for_site_group(
            curr_site_group_dict, cumulative_raw_per_sample_df,
            site_limited_per_sample_df, site_output_dir,
            curr_report_name, logger, config_fp=config_fp)

    return curr_per_sample_viral_loads_df, curr_cum_per_sample_base_df


def _extract_site_group_samples_only(per_sample_df, site_group):
    in_site_group_mask = per_sample_df[SITE_GROUP_KEY] == site_group
    site_group_samples_df = per_sample_df[in_site_group_mask].copy()
    return site_group_samples_df


def _is_context_site_group(site_group_name, config):
    is_context_site_group = False
    site_group_dict = get_site_group_dict_by_name(site_group_name, config)
    if site_group_dict is not None:
        is_context_site_group = site_group_dict.get(DO_CONTEXT_QC_KEY, False)
    return is_context_site_group


# if __name__ == "__main__":
#     generate_ucsf_report()
