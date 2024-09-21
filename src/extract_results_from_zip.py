import fnmatch
import glob
import os
import re
import zipfile
from pathlib import Path
from src.shared_symbols import *
from src.utils import extract_config_dict

ZIP_PATTERN = "*.zip"


def extract_to_folder(config_fp, zip_source_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    config = extract_config_dict(config_fp)

    _extract_results_from_zips(config, zip_source_dir, output_dir)


def _extract_results_from_zips(config_dict, inputs_dir, output_dir=None):
    output_dir = output_dir if output_dir is not None else inputs_dir
    qpcr_fname_pattern = config_dict[INPUT_FILES_KEY][FILENAME_PATTERN_KEY]
    qpcr_results_fps = glob.glob(f"{inputs_dir}/{ZIP_PATTERN}")
    qpcr_results_fps.sort(key=os.path.getmtime, reverse=False)
    for curr_qpcr_results_fp in qpcr_results_fps:
        print(f"INFO: extracting results file from {curr_qpcr_results_fp}")

        _extract_result_from_zip(
            curr_qpcr_results_fp, qpcr_fname_pattern, output_dir)


def _extract_result_from_zip(zip_fp, fname_pattern, output_dir):
    zip_file = zipfile.ZipFile(zip_fp)
    zip_filenames = zip_file.namelist()

    regex_pattern = fnmatch.translate(fname_pattern)
    raw_matching_fnames = [x for x in zip_filenames if
                           re.search(regex_pattern, x) is not None]
    matching_fnames = [x for x in raw_matching_fnames if
                       not x.startswith("__MACOSX") and "~" not in x]
    if len(matching_fnames) != 1:
        raise ValueError(f"Expected exactly one filename matching pattern "
                         f"'{fname_pattern}' in '{zip_fp}'")
    results_fname = matching_fnames[0]

    # extract a specific file from the zip container
    data = zip_file.read(results_fname)
    myfile_path = Path(output_dir) / Path(results_fname).name
    myfile_path.write_bytes(data)
    zip_file.close()

    return myfile_path
