import os
import pathlib
import shutil
from sys import argv
from src.shared_symbols import *
from src.utils import set_up_logger, extract_config_dict, \
    get_input_and_output_dirs
from src.extract_results_from_zip import extract_to_folder
from src.capture_qpcr_data import get_sorted_qpcr_results_fps


def intake_new_runs():
    """
    # example usage from command line:
     intake_new_runs \
       "/Users/abirmingham/Downloads" \
       "/Users/abirmingham/Work/Projects/covid_wastewater/qpcr/data/production" \
       "230609_UCSF_WW_H3K04ZF3_Results_20230609 135729.csv"
    """

    source_dir = argv[1]
    parent_output_dir = argv[2]
    datafile_fnames = argv[3:]

    logger, _ = set_up_logger()

    relocate_run_files(
        datafile_fnames, None, source_dir, parent_output_dir)


def relocate_run_files(
        datafile_fnames, config_fp, source_dir, parent_output_dir):
    for curr_datafile_fname in datafile_fnames:
        input_is_zip = pathlib.Path(curr_datafile_fname).suffix == ".zip"
        curr_full_run_name = pathlib.Path(curr_datafile_fname).stem

        curr_name_pieces = curr_full_run_name.split("_Results")
        curr_run_name = curr_name_pieces[0]

        source_fp = os.path.join(source_dir, curr_datafile_fname)
        run_dir = os.path.join(parent_output_dir, curr_run_name)
        new_fp = os.path.join(run_dir, curr_datafile_fname)

        if not os.path.exists(source_fp):
            raise ValueError(f"No file found named {source_fp}")

        # NOT using exist_ok=True ... if we already have a dir with this run
        # name, something is badly wrong
        os.makedirs(run_dir)

        # move the file
        os.rename(source_fp, new_fp)

        if input_is_zip:
            # run extract_results_from_zip on that zip file
            extract_to_folder(config_fp, run_dir, run_dir)
        # endif we need to extract the zip
    # next zip to pull


def centralize_run_files():
    """
    # example usage from command line:
     centralize_run_files \
       "/Users/abirmingham/Work/Projects/covid_wastewater/qpcr/data/production" \
       "/Users/abirmingham/Work/Projects/covid_wastewater/qpcr/data/temp/all_prods"
    """

    parent_source_dir = argv[1]
    new_dir_name = argv[2]

    parent_path = pathlib.Path(parent_source_dir)
    run_dirs = [str(x) for x in parent_path.iterdir() if x.is_dir()]
    centralize_run_files_from_dirs([None, run_dirs, new_dir_name])


def centralize_run_files_from_dirs(args_list):
    config_fp = args_list[0]
    source_dirs = args_list[1]
    parent_output_dir = args_list[2]

    config = extract_config_dict(config_fp)
    to_copy_fps = []
    for curr_source_dir in source_dirs:
        curr_results_fps = get_sorted_qpcr_results_fps(curr_source_dir, config)
        to_copy_fps.extend(curr_results_fps)

        curr_exclude_fps = get_sorted_qpcr_results_fps(
            curr_source_dir, config, EXCLUDE_FNAME_PATTERN_KEY)
        to_copy_fps.extend(curr_exclude_fps)
    # next source dir

    if len(to_copy_fps) > 0:
        # NB: our "output_dir" is the "input_dir" for the rest of the workflow
        output_dir, _ = get_input_and_output_dirs(parent_output_dir)
        os.makedirs(output_dir)
        for curr_fp in to_copy_fps:
            shutil.copy2(curr_fp, output_dir)
