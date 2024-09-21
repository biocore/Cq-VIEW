import logging
import pandas
import pathlib
import pprint
import os
import urllib.request
import yaml
from datetime import datetime
from src.shared_symbols import *


def get_input_and_output_dirs(source_dir):
    inputs_dir = os.path.join(source_dir, DATA_DIR_NAME)
    output_dir = os.path.join(source_dir, OUTPUT_DIR_NAME)
    return inputs_dir, output_dir


def extract_config_dict(config_fp=None):
    if config_fp is None:
        parent_dir = get_parent_dir()
        config_fp = os.path.join(parent_dir, "config.yml")

    # read in config file
    config_dict = extract_yaml_dict(config_fp)
    return config_dict


def extract_yaml_dict(yaml_fp):
    with open(yaml_fp, "r") as f:
        yaml_dict = yaml.safe_load(f)
    return yaml_dict


def get_assay_dict_by_name(assay_name, config):
    for curr_assay_dict in config[INPUT_FILES_KEY][ASSAYS_KEY]:
        if curr_assay_dict[ASSAY_NAME_KEY] == assay_name:
            return curr_assay_dict

    # if we didn't find an assay dict with the input name, return None
    return None


def extract_target_and_control_names(assay_dict):
    target_names = []
    control_names = []

    assay_targets_list = assay_dict[TARGETS_KEY]
    for curr_target_dict in assay_targets_list:
        curr_target_name = curr_target_dict[TARGET_NAME_KEY]
        curr_is_control_target = curr_target_dict.get(
            IS_FECAL_CONTROL_KEY, False)
        if curr_is_control_target:
            control_names.append(curr_target_name)
        else:
            target_names.append(curr_target_name)
    # next target in assay

    if len(target_names) != 1:
        raise ValueError(f"Expected exactly one target name but found "
                         f"'{target_names}'")

    if len(control_names) > 1:
        raise ValueError(f"Expected no more than one control name but found "
                         f"'{control_names}'")

    return target_names[0], control_names[0]


def get_site_group_dict_by_name(site_group_name, config):
    for curr_site_group_dict in config[SITE_GROUPS_KEY]:
        if curr_site_group_dict[SITE_GROUP_KEY] == site_group_name:
            return curr_site_group_dict

    # if we didn't find an assay dict with the input name, return None
    return None


def generate_timestamp_str():
    return datetime.now().strftime('%Y-%m-%d_%H-%M-%S')


def set_up_logger(log_fp=None, config=None):
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger()
    filehandler = None
    if log_fp is not None:
        filehandler = create_file_log_handler(log_fp)
        logger.addHandler(filehandler)

    if config is not None:
        logger.info("Configuration values:")
        logger.info(pprint.pformat(config))

    return logger, filehandler


def create_file_log_handler(log_fp):
    f_handler = logging.FileHandler(log_fp)
    f_handler.setLevel(logging.INFO)
    # f_format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    f_format = logging.Formatter('%(levelname)s: %(message)s')
    f_handler.setFormatter(f_format)
    return f_handler


def sort_df_by_date(input_df, date_col_key):
    # sort by date.  Has to be cast from string to real date to ensure
    # e.g. 7/7 sorts before 7/18
    temp_col_key = 'cast_date'
    input_df[temp_col_key] = pandas.to_datetime(input_df[date_col_key])
    input_df.sort_values(by=[temp_col_key], inplace=True)
    input_df.drop([temp_col_key], axis=1, errors="ignore", inplace=True)
    return input_df


def download_urls_to_dir(urls_fp, output_dir):
    def _add_end_backslash(a_str):
        if not a_str.endswith("/"):
            a_str = a_str + "/"
        return a_str

    output_dir = _add_end_backslash(output_dir)
    output_path = pathlib.Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # download all the urls in the input file, if there is one
    urls = []
    if urls_fp:
        with open(urls_fp) as urls_fh:
            # make sure to remove linebreaks ...
            urls = [x.strip() for x in urls_fh.readlines()]

    for curr_url in urls:
        curr_filename = pathlib.Path(curr_url).name
        curr_local_fp = output_path / curr_filename
        urllib.request.urlretrieve(curr_url, curr_local_fp)


def get_ref_dir():
    parent_dir = get_parent_dir()
    ref_dir = os.path.abspath(os.path.join(parent_dir, "reference_files"))
    return ref_dir


def get_parent_dir():
    curr_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.join(curr_dir, os.pardir)
    return parent_dir


def download_inputs_urls(output_dir, urls_fp=None):
    if urls_fp is None:
        ref_dir = get_ref_dir()
        urls_fp = os.path.join(ref_dir, "inputs_urls.txt")

    download_urls_to_dir(urls_fp, output_dir)
