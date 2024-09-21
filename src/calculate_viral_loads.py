import logging
import numpy
import pandas
import pprint
from datetime import timedelta

from src.shared_symbols import *
from src.utils import extract_config_dict, extract_target_and_control_names, \
    get_assay_dict_by_name
from src.capture_qpcr_data import extract_records_measured_for_target, \
    convert_measured_but_not_found_to_nan, get_cq_conf_key_for_target
from src.capture_and_plot_ladders import extract_ladders_for_target, \
    get_ladder_vals_from_config, get_mixed_model_intercepts_per_run_name, \
    get_single_curve_dicts_by_run
from src.output_formatted_data import output_campus_df, \
    output_standard_df_from_qpcr

RATIO_TIMES_TARGET_CONC_KEY = "old_norm_conc_is_ratio_times_target_conc"
NUM_REPS_USED_COL_NAME = 'num_replicates_used'
NUM_REPS_COL_NAME = 'num_replicates'

logging.basicConfig(level=logging.INFO)
# Must be defined here to be global
CONFIG = {}


def temp_ucsf_processing(sample_df, site_output_dir, run_name_w_site):
    pmmov_col = "PMMoV_N1_conc_per_L_from_legacy_std_curve"
    promega_col = "Promega_N1_conc_per_L_from_legacy_std_curve"

    # Go to sample-base df AFTER calculations; average the N1 and PMMoV
    # copies/L
    # NB: pivot_table, by default, averages the values for duplicate
    # indices, which is the desired behavior here.
    sample_base_mean_viral_loads_df = sample_df.pivot_table(
        index=SAMPLE_BASE_KEY,
        values=[pmmov_col, promega_col],
        dropna=False)
    sample_base_mean_viral_loads_df.reset_index(inplace=True)

    sample_base_mean_viral_loads_df[TARGET_TO_CONTROL_RATIO_KEY] = \
        sample_base_mean_viral_loads_df[
            promega_col] / sample_base_mean_viral_loads_df[pmmov_col] / 100

    info_df = _add_num_replicates_to_sample_base_df(
        sample_df, sample_base_mean_viral_loads_df, NUM_REPS_COL_NAME)
    info_df.sort_values(by=[SAMPLE_BASE_KEY], inplace=True)
    output_df = info_df[[SAMPLE_BASE_KEY, INTERNAL_COLLECTION_DATE_KEY,
                        TARGET_TO_CONTROL_RATIO_KEY]]

    output_df.to_csv(
        f"{site_output_dir}/{run_name_w_site}_delivery.csv", index=False)


def calc_viral_loads_for_site_group(
        site_group_dict, cumulative_raw_per_sample_df,
        site_limited_per_sample_df, output_dir, report_name, logger,
        config_fp=None):
    config_dict = extract_config_dict(config_fp)
    CONFIG.update(config_dict)

    do_context_qc = site_group_dict[DO_CONTEXT_QC_KEY]
    output_assay_dict = CONFIG[OUTPUT_FILES_KEY][OUTPUT_ASSAY_KEY]
    assay_name = output_assay_dict[ASSAY_NAME_KEY]

    assay_computed_curve_dicts = _get_computed_curve_dicts(
        cumulative_raw_per_sample_df, assay_name)

    # calculate per-sample viral loads & qc for that assay
    per_sample_viral_loads_df = _calc_per_sample_viral_loads_for_assay(
        site_limited_per_sample_df, assay_name, assay_computed_curve_dicts,
        do_context_qc)
    per_sample_viral_loads_df[OUTPUT_ID_KEY] = report_name
    per_sample_viral_loads_df[RESULT_ID_KEY] = \
        per_sample_viral_loads_df[INTERNAL_SAMPLE_AND_WELL_KEY] + "_" + \
        per_sample_viral_loads_df[OUTPUT_ID_KEY]
    per_sample_viral_loads_df.to_csv(
        f"{output_dir}/{report_name}_sample_data.csv", index=False)

    assay_dict = get_assay_dict_by_name(assay_name, CONFIG)
    target_name, _ = extract_target_and_control_names(assay_dict)
    output_campus_df(per_sample_viral_loads_df, target_name, output_dir,
                     report_name, CONFIG)
    output_standard_df_from_qpcr(
        per_sample_viral_loads_df, assay_name, target_name,
        output_dir, report_name, CONFIG)

    # calculate per-sample-base viral loads & qc for that assay
    per_sample_base_w_viral_loads_df = \
        _calc_per_sample_base_viral_loads_for_assay(
            per_sample_viral_loads_df, target_name)
    per_sample_base_w_viral_loads_df[OUTPUT_ID_KEY] = report_name

    sample_base_dfs = []
    if do_context_qc:
        site_prefixes = site_group_dict[SITE_PREFIXES_KEY]
        for curr_site_prefix in site_prefixes:
            logger.info(f"Sample_base context QC calculated for "
                        f"{curr_site_prefix}")
            curr_per_sample_base_w_viral_loads_w_qc_df = \
                _do_sample_base_w_context_qc_for_location(
                    per_sample_base_w_viral_loads_df,
                    curr_site_prefix, INTERNAL_NORM_CONC_KEY)

            curr_per_sample_base_w_viral_loads_w_qc_df.to_csv(
                f"{output_dir}/{report_name}_"
                f"sample_base_data_{curr_site_prefix}.csv",
                index=False)
            sample_base_dfs.append(curr_per_sample_base_w_viral_loads_w_qc_df)
        # next location

        cum_per_sample_base_df = pandas.concat(
            sample_base_dfs, ignore_index=True)
    else:
        per_sample_base_w_viral_loads_df.to_csv(
            f"{output_dir}/{report_name}_sample_base_data.csv",
            index=False)
        cum_per_sample_base_df = per_sample_base_w_viral_loads_df
    # endif

    return per_sample_viral_loads_df, cum_per_sample_base_df


def _get_computed_curve_dicts(cumulative_raw_per_sample_df, assay_name):
    assay_computed_curve_dicts = None
    curve_model_val = CONFIG[ANALYSIS_SETTINGS_KEY][CURVE_MODEL_KEY]
    if curve_model_val != POOLED_CURVE_VAL:
        assay_single_curves, assay_mixed_curves = \
            _get_single_and_mixed_curve_info_per_target(
                cumulative_raw_per_sample_df, assay_name)

        logger = logging.getLogger()
        if curve_model_val == MIXED_CURVE_VAL:
            logger.info(f"Mixed model curves for {assay_name} assay:")
            logger.info(pprint.pformat(assay_mixed_curves))
            assay_computed_curve_dicts = assay_mixed_curves
        elif curve_model_val == SINGLE_CURVE_VAL:
            logger.info(f"Single model curves for {assay_name} assay")
            logger.info(pprint.pformat(assay_single_curves))
            assay_computed_curve_dicts = assay_single_curves
        else:
            raise ValueError(f"Unrecognized model value '{curve_model_val}'")
        # end if curve model val is mixed/single
    # endif curve model is not pooled

    return assay_computed_curve_dicts


def _get_single_and_mixed_curve_info_per_target(cum_raw_df, assay_name):
    _, _, intercept_ladder_level = get_ladder_vals_from_config(CONFIG)

    intercept_info_per_target = {}
    single_curve_info_per_target = {}
    mixed_curve_info_per_target = {}
    assay_dict = get_assay_dict_by_name(assay_name, CONFIG)
    assay_targets = assay_dict[TARGETS_KEY]
    for curr_target_dict in assay_targets:
        curr_target_name = curr_target_dict[TARGET_NAME_KEY]

        # get all the ladder values for that target across days
        curr_ladder_df = extract_ladders_for_target(
            cum_raw_df, curr_target_dict, CONFIG)

        intercept_by_run_name = get_mixed_model_intercepts_per_run_name(
            curr_ladder_df, curr_target_name, CONFIG)
        intercept_info_per_target[curr_target_name] = intercept_by_run_name

        single_curve_dict_per_run, mixed_curve_dict_per_run = \
            get_single_curve_dicts_by_run(curr_ladder_df, curr_target_name,
                                          CONFIG, intercept_by_run_name)
        single_curve_info_per_target[curr_target_name] = \
            single_curve_dict_per_run
        mixed_curve_info_per_target[curr_target_name] = \
            mixed_curve_dict_per_run
    # next target

    return single_curve_info_per_target, mixed_curve_info_per_target


def _calc_per_sample_viral_loads_for_assay(
        site_limited_per_sample_df, assay_name, assay_computed_curve_dicts,
        do_context_qc):

    assay_and_site_limited_per_sample_df, target_name, control_name = \
        _extract_assay_samples_df_and_targets(
            site_limited_per_sample_df, assay_name)

    assay_dict = get_assay_dict_by_name(assay_name, CONFIG)
    assay_and_site_limited_per_sample_df[IGNORE_KEY] = False
    # TODO: Take out "do_context_qc as missing_target_fails_qc" hack.
    #  It is the case that the one site group that does
    #  context QC (san diego public wastewater) is ALSO the one site group
    #  where any missing target measurements should be treated as QC failures,
    #  but that is not guaranteed to be the case for every future
    #  site group, so we should specify these two setting separately rather
    #  than using do_context_qc *as* the missing_target_fails_qc value.
    assay_and_site_limited_per_sample_df = \
        _do_per_sample_qc(assay_and_site_limited_per_sample_df, assay_dict,
                          do_context_qc)

    # generate per-sample viral loads
    assay_and_site_limited_per_sample_w_viral_loads_df, _, fecal_control_curve_conc_key = \
        _calc_viral_loads_from_dilution_curves(
            assay_and_site_limited_per_sample_df, target_name, control_name,
            assay_computed_curve_dicts)
    assay_and_site_limited_per_sample_w_viral_loads_df.sort_values(
        by=[LOCATION_KEY, INTERNAL_COLLECTION_DATE_KEY], inplace=True)

    norm_model_val = CONFIG[ANALYSIS_SETTINGS_KEY][NORM_MODEL_KEY]
    if norm_model_val == NORM_SPIKEIN_VAL:
        per_sample_normed_df = \
            assay_and_site_limited_per_sample_w_viral_loads_df
        per_sample_normed_df.rename(
            columns={RATIO_TIMES_TARGET_CONC_KEY: INTERNAL_NORM_CONC_KEY},
            inplace=True)
    elif norm_model_val == NORM_ROLLING_VAL:
        per_sample_normed_df, rolling_avg_key = _calc_rolling_avg(
            fecal_control_curve_conc_key,
            assay_and_site_limited_per_sample_w_viral_loads_df)
        per_sample_normed_df[INTERNAL_NORM_CONC_KEY] = \
            per_sample_normed_df.apply(
                lambda row: _calc_single_rolling_avg_norm_conc_by_location(
                    row[TARGET_TO_CONTROL_RATIO_KEY], row[rolling_avg_key],
                    row[LOCATION_KEY]),
            axis=1)
    elif norm_model_val == NORM_FIXED_VAL:
        per_sample_normed_df = \
            assay_and_site_limited_per_sample_w_viral_loads_df
        per_sample_normed_df[INTERNAL_NORM_CONC_KEY] = \
            per_sample_normed_df.apply(
                lambda row: _calc_single_fixed_norm_conc_by_location(
                    row[TARGET_TO_CONTROL_RATIO_KEY], row[LOCATION_KEY]),
                axis=1)
    else:
        raise ValueError(f"Unrecognized normalization model "
                         f"'{norm_model_val}'")

    # do qc that is appropriate to this assay and target
    # (done here not earlier bc may need to use per-sample viral loads)
    per_sample_normed_df = \
        _do_per_sample_context_and_base_qc(
            per_sample_normed_df, assay_name, INTERNAL_NORM_CONC_KEY,
            do_context_qc)

    return per_sample_normed_df


def _calc_single_fixed_norm_conc_by_location(
        target_to_control_ratio, location_prefix):

    # look up the fixed norm factor for that prefix in the analysis settings;
    # if there isn't one, use default
    fixed_norm_factors = CONFIG[ANALYSIS_SETTINGS_KEY][NORM_FIXED_FACTORS]
    norm_factor_for_location = fixed_norm_factors.get(
        location_prefix, fixed_norm_factors["default"])

    return target_to_control_ratio * norm_factor_for_location


def _calc_single_rolling_avg_norm_conc_by_location(
        target_to_control_ratio, rolling_avg, location_prefix):

    # look up the divisor for that prefix in the analysis settings;
    # if there isn't one, use NaN
    divisors = CONFIG[ANALYSIS_SETTINGS_KEY][NORM_ROLLING_DIVISORS]
    divisor_for_location = divisors.get(location_prefix, numpy.NaN)

    return target_to_control_ratio * rolling_avg / divisor_for_location


def _extract_assay_samples_df_and_targets(per_samples_df, assay_name):
    assay_per_sample_df = None
    assay_dict = get_assay_dict_by_name(assay_name, CONFIG)
    target_name, control_name = extract_target_and_control_names(assay_dict)
    logger = logging.getLogger()

    if control_name is not None and \
            control_name not in per_samples_df.columns:
        logger.info(f"Control '{control_name}' not in qpcr data")
        return assay_per_sample_df, target_name, control_name

    if target_name not in per_samples_df.columns:
        logger.info(f"Target '{target_name}' not in qpcr data")
        return assay_per_sample_df, target_name, control_name

    # get only records from the input df that are not NaN for target
    assay_per_sample_df = extract_records_measured_for_target(
        per_samples_df, target_name)
    if control_name is not None:
        lacks_control_measurement_mask = \
            assay_per_sample_df[control_name].isna()
        if lacks_control_measurement_mask.any():
            raise ValueError(f"Some records were measured for target "
                             f"'{target_name}' but not for control "
                             f"'{control_name}'")

    # convert internal_not_measured value into nan
    assay_per_sample_df = convert_measured_but_not_found_to_nan(
        assay_per_sample_df)
    return assay_per_sample_df, target_name, control_name


def _do_per_sample_qc(per_sample_df, assay_dict, missing_target_fails_qc):
    per_sample_df[FAILS_SAMPLE_QC_KEY] = False
    curr_targets_dict = assay_dict[TARGETS_KEY]
    for curr_dict in curr_targets_dict:
        if TARGET_LADDER_PREFIX_KEY in curr_dict:
            curr_is_target = not curr_dict.get(IS_FECAL_CONTROL_KEY, False)
            if curr_is_target:
                curr_has_na_target_mask = \
                    _get_target_is_na_mask(per_sample_df, curr_dict)
                per_sample_df[TARGET_PRESENT_KEY] = ~curr_has_na_target_mask
            # endif this is the assay target

            per_sample_df = _do_cq_conf_check(per_sample_df, curr_dict)

            # if this is the target gene, and if this data should NOT fail
            # in the case of the target gene not being found, then skip
            # the QC to check if this gene is found at all.
            # (For example, for samples from individual buildings, it should
            # NOT be a failure if the target virus is not found, because
            # it is realistic that there's just no one in that building who is
            # infected.  However, for samples from wastewater treatment plants,
            # it is not realistic that *no one* in the catchment area is
            # infected, so in that case we mark any cases where the target
            # is not found at all as QC failures.)
            if curr_is_target and not missing_target_fails_qc:
                continue
            per_sample_df = _do_no_assay_nas_check(per_sample_df, curr_dict)

        # endif this target has ladder prefixes
    # next target in this assay

    per_sample_df = _do_no_exclude_note_check(per_sample_df)
    return per_sample_df


def _do_no_exclude_note_check(per_sample_df):
    has_exclude_note_mask = per_sample_df[EXCLUDE_NOTE_KEY].notna()
    curr_sample_w_check_true = per_sample_df.loc[
            has_exclude_note_mask,
            INTERNAL_SAMPLE_AND_WELL_KEY].tolist()

    per_sample_df = _do_specified_qc_check(
        per_sample_df, curr_sample_w_check_true, EXCLUDE_NOTE_KEY,
        "present", "Has an exclude note", FAILS_SAMPLE_QC_KEY)

    return per_sample_df


def _do_no_assay_nas_check(per_sample_df, target_dict):
    curr_has_na_target_mask = \
        _get_target_is_na_mask(per_sample_df, target_dict)
    curr_sample_w_check_true = per_sample_df.loc[
            curr_has_na_target_mask,
            INTERNAL_SAMPLE_AND_WELL_KEY].tolist()

    curr_target = target_dict[TARGET_NAME_KEY]
    per_sample_df = _do_specified_qc_check(
        per_sample_df, curr_sample_w_check_true, curr_target,
        "na", "Has NA", FAILS_SAMPLE_QC_KEY)

    return per_sample_df


def _get_target_is_na_mask(per_sample_df, target_dict):
    curr_target = target_dict[TARGET_NAME_KEY]
    curr_has_na_target_mask = per_sample_df[curr_target].isna()
    return curr_has_na_target_mask


def _do_cq_conf_check(per_sample_df, target_dict, cq_gte=None):
    if cq_gte is None:
        cq_gte = CONFIG[ANALYSIS_SETTINGS_KEY][MIN_ACCEPTABLE_CQ_CONF_KEY]

    curr_target_name = target_dict[TARGET_NAME_KEY]
    curr_cq_conf_key = get_cq_conf_key_for_target(
        curr_target_name, CONFIG)

    # Note the limitation that this method is only qc-ing cq_confidence values
    # *greater than* zero.  That's because when a target is not found
    # (Undetermined), its cq confidence is zero--and in no other situation,
    # it appears.  Since qc-ing presence/absence of the target is outside the
    # scope of this method, we only check the cq_confidence values
    # for being too low if they are greater than zero.
    curr_has_cq_gt_zero_mask = per_sample_df[curr_cq_conf_key] > 0
    curr_has_low_cq_mask = per_sample_df[curr_cq_conf_key] < cq_gte
    curr_combined_mask = curr_has_cq_gt_zero_mask & curr_has_low_cq_mask
    curr_sample_w_check_true = per_sample_df.loc[
        curr_combined_mask, INTERNAL_SAMPLE_AND_WELL_KEY].tolist()

    per_sample_df = _do_specified_qc_check(
        per_sample_df, curr_sample_w_check_true, curr_target_name,
        "cq_conf", f"Has cq confidence < {cq_gte}",
        FAILS_SAMPLE_QC_KEY)

    return per_sample_df


def _do_specified_qc_check(a_df, items_w_check_true_list, target_name,
                           check_name, check_note, qc_level_key):
    a_df, target_check_key = _set_up_check_column(
        a_df, target_name, check_name, None)

    if qc_level_key in [FAILS_SAMPLE_BASE_QC_KEY,
                        FAILS_SAMPLE_BASE_CONTEXT_QC_KEY]:
        item_key = SAMPLE_BASE_KEY
    elif qc_level_key in [FAILS_SAMPLE_QC_KEY, FAILS_SAMPLE_CONTEXT_QC_KEY]:
        item_key = INTERNAL_SAMPLE_AND_WELL_KEY
    else:
        raise ValueError(f"No item key known for qc key '{qc_level_key}'")

    has_check_true_sample_mask = a_df[item_key].isin(
        items_w_check_true_list)
    a_df[target_check_key] = None
    a_df.loc[has_check_true_sample_mask, target_check_key] = \
        f"{check_note} for {target_name}"
    a_df[qc_level_key] = \
        a_df[qc_level_key] | has_check_true_sample_mask
    a_df[IGNORE_KEY] = \
        a_df[IGNORE_KEY] | has_check_true_sample_mask

    return a_df


def _set_up_check_column(a_df, target_name, check_name, default_val=None):
    target_check_key = f"{target_name}_{check_name}_check"
    a_df[target_check_key] = default_val
    return a_df, target_check_key


def _calc_viral_loads_from_dilution_curves(
        working_df, target_name, fecal_control_name,
        assay_computed_curve_dicts):

    target_curves_per_run = _get_target_curve_dicts_per_run(
        assay_computed_curve_dicts, target_name, working_df)
    working_df, target_curve_conc_key = _calc_concentrations_from_curve(
        target_name, target_curves_per_run, working_df)

    target_curves_per_run = _get_target_curve_dicts_per_run(
        assay_computed_curve_dicts, fecal_control_name, working_df)
    working_df, fecal_control_curve_conc_key = _calc_concentrations_from_curve(
        fecal_control_name, target_curves_per_run, working_df)

    working_df[TARGET_TO_CONTROL_RATIO_KEY] = \
        working_df[target_curve_conc_key]/working_df[
            fecal_control_curve_conc_key]

    target_spikein_curve_dict = \
        CONFIG[SPIKEIN_CURVE_KEY].get(target_name, None)
    if target_spikein_curve_dict is None:
        logger = logging.getLogger()
        logger.info(f"No spike-in for target {target_name}")
        # pivot_values = TARGET_TO_CONTROL_RATIO_KEY
    else:
        # pivot_values = INTERNAL_NORM_CONC_KEY
        spikein_curves_per_run = _get_target_curve_dicts_per_run(
            None, target_name, working_df, get_spikein=True)
        working_df, target_spikein_conc_key = _calc_concentrations_from_curve(
            target_name, spikein_curves_per_run, working_df)

        working_df[RATIO_TIMES_TARGET_CONC_KEY] = \
            working_df[TARGET_TO_CONTROL_RATIO_KEY] * \
            working_df[target_spikein_conc_key]

    # TODO: if removing spikein calc for good, modify return values
    return working_df, None, fecal_control_curve_conc_key


def _get_target_curve_dicts_per_run(
        curves_per_target_per_run, target_name, working_df, get_spikein=False):
    unique_run_names = working_df[RUN_NAME_KEY].drop_duplicates()
    if curves_per_target_per_run is None:
        curve_key = SPIKEIN_CURVE_KEY if get_spikein else STD_CURVE_KEY
        target_curves_per_run = \
            {x: CONFIG[curve_key][target_name] for x in unique_run_names}
    else:
        target_curves_per_run = curves_per_target_per_run[target_name]
    return target_curves_per_run


def _calc_concentrations_from_curve(
        target_name, curve_dict_by_run_name, working_df):
    first_run_name = list(curve_dict_by_run_name.keys())[0]
    first_run_curve = curve_dict_by_run_name[first_run_name]
    curve_name = first_run_curve[CURVE_NAME_KEY]
    new_conc_key = f"{target_name}{CONC_SUFFIX}{curve_name}"

    working_df[new_conc_key] = working_df.apply(
        lambda row: _calc_single_conc_from_curve(
            curve_dict_by_run_name, row[target_name], row[RUN_NAME_KEY]),
        axis=1)
    return working_df, new_conc_key


def _calc_single_conc_from_curve(
        target_curve_dict_by_run_name, raw_target_cq, run_name):

    target_curve_dict = target_curve_dict_by_run_name[run_name]
    target_conc = _calc_single_raw_conc_per_l_from_std_curve(
        raw_cq=raw_target_cq, intercept=target_curve_dict[INTERCEPT_KEY],
        slope=target_curve_dict[SLOPE_KEY],
        reaction_vol_per_l=target_curve_dict[REACTION_VOL_PER_L_KEY],
        is_natural_log=target_curve_dict[USE_LN_KEY])
    return target_conc


# NB: * as first arg means no position args allowed--all args must be named
def _calc_single_raw_conc_per_l_from_std_curve(*, raw_cq, intercept, slope,
                                               reaction_vol_per_l,
                                               is_natural_log=True):
    # NB: "raw" here means un-normalized.
    # A standard curve has the form
    # raw target Cq =
    # slope * [ln or log10] (target gene copies per some volume) + intercept
    # Since we have raw target Cq and want to calculate target gene copies/L,
    # (the target concentration), this equation must be reformatted to:
    # target gene copies/L =
    # [e or 10] ^ ((raw target Cq - intercept)/slope) * volume_multiplier
    # where reaction_vol_per_l is what we have to multiply the volume used for
    # the standard curve by to get to liters (for example, if the standard
    # curve was done in microliters, we have to multiply that by 1 million to
    # get to liters.)
    # Note that the target gene copies/L is also "raw"--unnormalized.

    exponent_term = (raw_cq - intercept)/slope
    if is_natural_log:
        raw_target_conc_per_reaction = numpy.exp(exponent_term)
    else:  # is "common" log--ie., base 10
        raw_target_conc_per_reaction = 10 ** exponent_term

    raw_target_conc_per_l = raw_target_conc_per_reaction * reaction_vol_per_l
    return raw_target_conc_per_l


def _calc_rolling_avg(col_name, working_df):
    rolling_avg_key = f"{col_name}_rolling_avg_by_location"
    working_df[rolling_avg_key] = working_df.apply(
        lambda row: _calc_rolling_avg_from_date(
            working_df, col_name, row[LOCATION_KEY],
            row[INTERNAL_COLLECTION_DATE_KEY]),
        axis=1)

    return working_df, rolling_avg_key


def _calc_rolling_avg_from_date(
        per_sample_df, col_name, a_location, a_date):

    rolling_avg_num_days = CONFIG[ANALYSIS_SETTINGS_KEY][ROLLING_AVG_DAYS_KEY]

    # make a mask:
    lte_curr_date_mask = per_sample_df[INTERNAL_COLLECTION_DATE_KEY] <= a_date
    gt_threshold_date_mask = \
        per_sample_df[INTERNAL_COLLECTION_DATE_KEY] > a_date - \
        timedelta(days=rolling_avg_num_days)
    at_desired_location_mask = per_sample_df[LOCATION_KEY] == a_location
    passes_indiv_qc_mask = per_sample_df[FAILS_SAMPLE_QC_KEY] == False # noqa E712
    relevant_rows_mask = \
        lte_curr_date_mask & gt_threshold_date_mask & \
        at_desired_location_mask & passes_indiv_qc_mask

    # get all the values in col_name for rows that match the mask
    unique_relevant_dates = pandas.unique(
        per_sample_df.loc[relevant_rows_mask, INTERNAL_COLLECTION_DATE_KEY])
    rolling_avg_min_num_days = \
        CONFIG[ANALYSIS_SETTINGS_KEY][ROLLING_AVG_MIN_DAYS_KEY]
    if len(unique_relevant_dates) < rolling_avg_min_num_days:
        return numpy.NaN

    # take avg of all the target_name concs for rows that match the mask
    relevant_col_vals = per_sample_df.loc[relevant_rows_mask, col_name]
    avg_of_relevant_vals = relevant_col_vals.mean()
    return avg_of_relevant_vals


def _do_per_sample_context_and_base_qc(
        per_sample_df, assay_name, per_sample_viral_load_key, do_context_qc):

    # TODO: should I do the 2-replicate test on the raw data too?

    per_sample_df[FAILS_SAMPLE_CONTEXT_QC_KEY] = False
    if do_context_qc:
        rep_fold_change_threshold = \
            CONFIG[ANALYSIS_SETTINGS_KEY][REPLICATE_FOLD_CHANGE_THRESHOLD_KEY]
        per_sample_df = _do_too_large_ratio_for_multiple_replicates_check(
            per_sample_df, per_sample_viral_load_key,
            rep_fold_change_threshold)

    per_sample_df[FAILS_SAMPLE_BASE_QC_KEY] = False
    min_num_passing_reps = CONFIG[ANALYSIS_SETTINGS_KEY][MIN_PASSING_REPS]
    passes_per_sample_qc_mask = per_sample_df[IGNORE_KEY] == False # noqa E712
    per_sample_df = _do_min_replicates_check(
        per_sample_df, assay_name, min_num_passing_reps,
        "passing", passes_per_sample_qc_mask)

    return per_sample_df


def _do_too_large_ratio_for_multiple_replicates_check(
        per_sample_df, target_name, fold_change_threshold):

    all_failing_indexes = []
    passes_per_sample_qc_mask = per_sample_df[IGNORE_KEY] == False # noqa E712
    passing_per_sample_df = per_sample_df.loc[passes_per_sample_qc_mask, :]

    sample_base_grouping = passing_per_sample_df.groupby(SAMPLE_BASE_KEY)
    for curr_sample_base_key, curr_sample_base_df in sample_base_grouping:
        curr_failing_indexes = \
            _get_indexes_of_records_failing_fold_change_check(
                curr_sample_base_df, target_name, fold_change_threshold)
        all_failing_indexes.extend(curr_failing_indexes)
    # next group (i.e., sample_base)

    samples_failing_check = \
        per_sample_df.loc[all_failing_indexes,
                          INTERNAL_SAMPLE_AND_WELL_KEY].tolist()
    per_sample_df = _do_specified_qc_check(
        per_sample_df, samples_failing_check, target_name,
        "rep_ratio", f"Has {target_name} >= {fold_change_threshold} "
                     f"of most distant accepted replicate",
        FAILS_SAMPLE_CONTEXT_QC_KEY)
    return per_sample_df


def _get_indexes_of_records_failing_fold_change_check(
        sample_base_df, target_name, fold_change_threshold):

    sample_base_vals = sample_base_df[target_name].copy()
    sample_base_median = sample_base_vals.median()
    diff_from_base_median = sample_base_vals - sample_base_median
    abs_diff_from_base_median = diff_from_base_median.abs()
    abs_diff_from_base_median.sort_values(inplace=True)

    curr_passing_indexes = []
    for curr_abs_idx, curr_abs_diff in abs_diff_from_base_median.items():
        curr_expanded_indexes = curr_passing_indexes.copy()
        curr_expanded_indexes.append(curr_abs_idx)
        if len(curr_expanded_indexes) >= 2:
            vals_to_check = \
                sample_base_df.loc[curr_expanded_indexes, target_name]
            curr_fold_change = _get_fold_change_ratio(vals_to_check)
            if curr_fold_change > fold_change_threshold:
                break
            # endif fold_change > threshold
        # endif len(curr_expanded_indexes) >=2

        curr_passing_indexes = curr_expanded_indexes.copy()
    # next curr_abs

    # if this sample has two or more replicates, but there is only one
    # replicate listed as passing the fc check, that means really nothing
    # passed the fc check
    if len(abs_diff_from_base_median) >= 2:
        if len(curr_passing_indexes) == 1:
            curr_passing_indexes = []

    failing_indexes = \
        list(set(abs_diff_from_base_median.index) - set(curr_passing_indexes))
    return failing_indexes


def _get_fold_change_ratio(vals_to_check):
    vals_max = vals_to_check.max()
    vals_min = vals_to_check.min()
    vals_ratio = vals_max / vals_min
    return vals_ratio


def _do_min_replicates_check(
        per_sample_df, assay_name, min_num_reps,
        addtl_condition_str="", addtl_condition_mask=None):

    REP_COUNT_KEY = "rep_count"
    working_df = per_sample_df
    if addtl_condition_mask is not None:
        working_df = per_sample_df.loc[addtl_condition_mask, :]

    per_sample_base_count = working_df[SAMPLE_BASE_KEY].value_counts(
        dropna=False)
    qc_per_sample_base_df = pandas.DataFrame(
        {REP_COUNT_KEY: per_sample_base_count})
    qc_per_sample_base_df.reset_index(inplace=True)
    qc_per_sample_base_df.rename(
        columns={"index": SAMPLE_BASE_KEY}, inplace=True)

    not_min_replicates_mask = \
        qc_per_sample_base_df[REP_COUNT_KEY] < min_num_reps
    not_min_replicates_sample_base_list = qc_per_sample_base_df.loc[
            not_min_replicates_mask, SAMPLE_BASE_KEY].tolist()

    check_name = REP_COUNT_KEY
    if addtl_condition_str:
        check_name = f"{addtl_condition_str}_{REP_COUNT_KEY}"

    per_sample_df = _do_specified_qc_check(
        per_sample_df, not_min_replicates_sample_base_list,
        f"{assay_name}_sample_base", check_name,
        f"Does not have gte {min_num_reps} {addtl_condition_str} replicates",
        FAILS_SAMPLE_BASE_QC_KEY)

    return per_sample_df


def _calc_per_sample_base_viral_loads_for_assay(
        per_sample_w_rolling_avg_df, target_name):
    input_df = per_sample_w_rolling_avg_df.copy()
    per_sample_base_w_viral_loads_df = \
        _build_up_complete_sample_base_df(
            input_df,
            INTERNAL_NORM_CONC_KEY, target_name)

    return per_sample_base_w_viral_loads_df


def _build_up_complete_sample_base_df(
        per_sample_df, target_norm_conc_key, target_cq_key):
    # get items with ignore = false and target present = true
    # for these, calc viral load (and num replicates used)
    good_pos_sample_bases_df, good_pos_sample_bases = \
        _get_relevant_sample_base_df(
            per_sample_df, target_norm_conc_key, target_cq_key,
            ignore_val=False, target_present_val=True)
    good_pos_sample_bases_df[IGNORE_KEY] = False

    # get items with ignore = false and target present = false AND
    # sample base is not in previous set;
    # for these, set viral load to -1 (measured but not found)
    good_neg_sample_bases_df, good_neg_sample_bases = \
        _get_relevant_sample_base_df(
            per_sample_df, target_norm_conc_key, target_cq_key,
            ignore_val=False, target_present_val=False,
            excluded_sample_bases=good_pos_sample_bases)
    good_neg_sample_bases_df[IGNORE_KEY] = False

    # get items with ignore = true AND sample base not in either of
    # previous sets then set viral load to NA for all
    good_sample_bases = good_pos_sample_bases + good_neg_sample_bases
    bad_sample_bases_df, bad_sample_bases = _get_relevant_sample_base_df(
        per_sample_df, target_norm_conc_key, target_cq_key,
        ignore_val=True, target_present_val=None,
        excluded_sample_bases=good_sample_bases)
    bad_sample_bases_df[IGNORE_KEY] = True

    known_sample_bases = good_sample_bases + bad_sample_bases
    df_sample_bases = pandas.unique(per_sample_df.loc[:, SAMPLE_BASE_KEY])
    if set(known_sample_bases) != set(df_sample_bases):
        # TODO: if keep this, give it a better message
        raise ValueError("Some base(s) missed when building complete df")

    complete_per_sample_base_df = pandas.concat(
        [good_pos_sample_bases_df, good_neg_sample_bases_df,
         bad_sample_bases_df], ignore_index=True)

    complete_per_sample_base_df[UNRELIABLE_FLAG] = \
        complete_per_sample_base_df[NUM_REPS_USED_COL_NAME] < 2
    complete_per_sample_base_df = _add_num_replicates_to_sample_base_df(
        per_sample_df, complete_per_sample_base_df, NUM_REPS_COL_NAME)

    return complete_per_sample_base_df


def _get_relevant_sample_base_df(
        per_sample_df, target_norm_conc_key, target_cq_key,
        ignore_val=False, target_present_val=None, excluded_sample_bases=None):

    if excluded_sample_bases is None:
        excluded_sample_bases = []

    ignore_val_mask = per_sample_df[IGNORE_KEY] == ignore_val
    not_excluded_sample_base_mask = \
        ~per_sample_df[SAMPLE_BASE_KEY].isin(excluded_sample_bases)
    combined_mask = ignore_val_mask & not_excluded_sample_base_mask
    if target_present_val is not None:
        target_present_val_mask = \
            per_sample_df[TARGET_PRESENT_KEY] == target_present_val
        combined_mask = combined_mask & target_present_val_mask

    relevant_per_sample_df = \
        per_sample_df.loc[combined_mask, :].copy()

    per_sample_base_w_viral_loads_df = \
        _get_per_sample_base_viral_load_df(
            relevant_per_sample_df, target_norm_conc_key, target_cq_key)

    if not ignore_val:
        if not target_present_val:
            per_sample_base_w_viral_loads_df[target_norm_conc_key] = \
                INTERNAL_MEASURED_BUT_NOT_FOUND_VALUE
            per_sample_base_w_viral_loads_df[target_cq_key] = \
                INTERNAL_MEASURED_BUT_NOT_FOUND_VALUE
        # endif not target present
    else:  # if ignore_val == True
        per_sample_base_w_viral_loads_df[target_norm_conc_key] = numpy.nan
        per_sample_base_w_viral_loads_df[target_cq_key] = numpy.nan
    # endif ignore_val == False/ == True

    relevant_sample_bases = \
        pandas.unique(per_sample_base_w_viral_loads_df.loc[:, SAMPLE_BASE_KEY])

    return per_sample_base_w_viral_loads_df, list(relevant_sample_bases)


def _get_per_sample_base_viral_load_df(
        sample_df, per_sample_viral_load_key, target_cq_key):

    # Go to sample-base df AFTER calculations; average the viral loads.
    # NB: pivot_table, by default, averages the values for duplicate
    # indices, which is the desired behavior here.
    sample_base_mean_viral_loads_df = sample_df.pivot_table(
        index=SAMPLE_BASE_KEY,
        values=[per_sample_viral_load_key, target_cq_key], dropna=False)
    sample_base_mean_viral_loads_df.reset_index(inplace=True)

    output_df = _add_num_replicates_to_sample_base_df(
        sample_df, sample_base_mean_viral_loads_df)
    output_df = output_df[
        [SAMPLE_BASE_KEY, INTERNAL_COLLECTION_DATE_KEY,
         per_sample_viral_load_key, target_cq_key,
         LOCATION_KEY, NUM_REPS_USED_COL_NAME]]

    return output_df


def _add_num_replicates_to_sample_base_df(
        sample_df, sample_base_df, reps_col_name=NUM_REPS_USED_COL_NAME):
    base_location_date_sets_df = _get_num_replicates_per_sample_base(
        sample_df, reps_col_name)
    merge_on = [SAMPLE_BASE_KEY]
    for curr_col in [INTERNAL_COLLECTION_DATE_KEY, LOCATION_KEY]:
        if curr_col in sample_base_df.columns:
            merge_on.append(curr_col)
    output_df = sample_base_df.merge(
        base_location_date_sets_df, on=merge_on)
    output_df.sort_values(by=[INTERNAL_COLLECTION_DATE_KEY], inplace=True)
    return output_df


def _get_num_replicates_per_sample_base(
        per_sample_df, reps_col_name):
    base_location_date_sets_df = \
        per_sample_df.groupby([SAMPLE_BASE_KEY,
                               LOCATION_KEY,
                               INTERNAL_COLLECTION_DATE_KEY],
                              dropna=False).size().reset_index()
    base_location_date_sets_df.rename(
        columns={0: reps_col_name}, inplace=True)
    base_location_date_sets_df = base_location_date_sets_df[
        [SAMPLE_BASE_KEY, INTERNAL_COLLECTION_DATE_KEY, LOCATION_KEY,
         reps_col_name]]
    return base_location_date_sets_df


def _do_sample_base_w_context_qc_for_location(
        per_sample_base_df, location_prefix, per_sample_viral_load_key):

    temp_mask = per_sample_base_df[LOCATION_KEY] == location_prefix
    temp_df = per_sample_base_df.loc[temp_mask, :].copy()

    per_sample_base_w_viral_loads_w_qc_df = \
        _do_sample_base_context_qc(temp_df, per_sample_viral_load_key)

    return per_sample_base_w_viral_loads_w_qc_df


def _do_sample_base_context_qc(per_sample_base_df, per_sample_viral_load_key):
    per_sample_base_df[FAILS_SAMPLE_BASE_CONTEXT_QC_KEY] = False
    per_sample_base_df = _do_too_different_from_neighbors_check(
        per_sample_base_df, per_sample_viral_load_key)
    return per_sample_base_df


def _do_too_different_from_neighbors_check(per_sample_base_df, col_name):
    context_fold_change_threshold = \
        CONFIG[ANALYSIS_SETTINGS_KEY][CONTEXT_FOLD_CHANGE_THRESHOLD_KEY]

    # sort by date ascending (note this is making a copy) and re-index
    sorted_per_sample_df = per_sample_base_df.sort_values(
        by=[INTERNAL_COLLECTION_DATE_KEY])
    sorted_per_sample_df = \
        sorted_per_sample_df[sorted_per_sample_df[INTERNAL_COLLECTION_DATE_KEY].notna()]
    sorted_per_sample_df = \
        sorted_per_sample_df[sorted_per_sample_df[INTERNAL_NORM_CONC_KEY].notna()]
    sorted_per_sample_df.reset_index(drop=True, inplace=True)

    # if not every date is unique, error out
    unique_dates = pandas.unique(
        sorted_per_sample_df[INTERNAL_COLLECTION_DATE_KEY])
    if len(unique_dates) != len(sorted_per_sample_df):
        raise ValueError("Expected all dates to be unique")

    failing_sample_bases = []
    context_window = CONFIG[ANALYSIS_SETTINGS_KEY][CONTEXT_DAYS_WINDOW_KEY]
    for curr_idx in range(0, len(sorted_per_sample_df)):
        curr_row_dict = _get_row_dict_by_index(sorted_per_sample_df, curr_idx)
        context_fold_changes = []
        for curr_day_diff in range(1, context_window+1):
            prev_row_dict = _get_row_dict_by_index(
                sorted_per_sample_df, curr_idx - curr_day_diff)
            next_row_dict = _get_row_dict_by_index(
                sorted_per_sample_df, curr_idx + curr_day_diff)
            prev_fold_change = _calculate_fold_change_between_rows(
                curr_row_dict, prev_row_dict, col_name)
            next_fold_change = _calculate_fold_change_between_rows(
                curr_row_dict, next_row_dict, col_name)
            context_fold_changes.append(prev_fold_change)
            context_fold_changes.append(next_fold_change)
        # next day diff

        w_in_threshold_fcs = [x for x in context_fold_changes
                              if x is not None and
                              x <= context_fold_change_threshold]

        if len(w_in_threshold_fcs) == 0:
            failing_sample_bases.append(curr_row_dict[SAMPLE_BASE_KEY])

    per_sample_base_df = _do_specified_qc_check(
        per_sample_base_df, failing_sample_bases, col_name, "context",
        f"Is not within {context_fold_change_threshold}-fold "
        f"of a context value", FAILS_SAMPLE_BASE_CONTEXT_QC_KEY)

    return per_sample_base_df


def _get_row_dict_by_index(per_sample_base_df, an_idx):
    if an_idx < 0 or an_idx > len(per_sample_base_df)-1:
        return None

    idx_dict = per_sample_base_df.iloc[an_idx].to_dict()
    return idx_dict


def _calculate_fold_change_between_rows(
        curr_row_dict, other_row_dict, col_name):
    if curr_row_dict is None or other_row_dict is None:
        return None

    # TODO: Put back?  Check not to return value if it is > X days away from
    #  the value we are checking the context for ...
    # curr_row_date = curr_row_dict[INTERNAL_COLLECTION_DATE_KEY]
    # other_row_date = other_row_dict[INTERNAL_COLLECTION_DATE_KEY]
    # context_min_days_threshold = \
    #     CONFIG[ANALYSIS_SETTINGS_KEY][CONTEXT_DAYS_THRESHOLD_KEY]
    # if abs(curr_row_date - other_row_date) >
    #       timedelta(days=context_min_days_threshold):
    #     return None

    curr_row_col_val = curr_row_dict[col_name]
    other_row_col_val = other_row_dict[col_name]
    fold_change = (max(curr_row_col_val, other_row_col_val) /
                   min(curr_row_col_val, other_row_col_val))

    return fold_change
