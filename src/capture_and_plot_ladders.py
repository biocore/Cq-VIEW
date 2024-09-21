import logging
import numpy
import os
import pandas
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure
from matplotlib.ticker import MultipleLocator
from src.shared_symbols import *
from src.utils import set_up_logger, extract_config_dict, \
    generate_timestamp_str, get_input_and_output_dirs
from src.capture_qpcr_data import get_cq_conf_key_for_target, \
    capture_qpcr_data, combine_all_datesets, \
    convert_measured_but_not_found_to_nan, get_unique_sorted_date_tuples

# config file keys
LADDERS_KEY = "ladders"
LADDER_BOTTOM_KEY = "ladder_bottom"
LADDER_STEPS_KEY = "ladder_steps"
REGRESSION_BOTTOM_STEP_KEY = "regression_bottom_step"
LADDER_STEP_KEY = "ladder_step"
REPLICATE_KEY = "replicate"
INTERCEPT_LADDER_STEP_KEY = "intercept_ladder_step"

# config key added from inputs
INCLUDE_MIXED_KEY = "include_mixed"

# sample name replicate separator string
REPLICATE_SEPARATOR = ".R"


def capture_and_save_ladders(
        source_dirs, include_mixed_model_on_plots, config_fp,
        timestamp_str=None):

    try:
        config = extract_config_dict(config_fp)
        if not timestamp_str:
            timestamp_str = generate_timestamp_str()
        config[INCLUDE_MIXED_KEY] = include_mixed_model_on_plots

        counter = 1
        ladder_fps = []
        plot_fps = []
        for source_dir in source_dirs:
            inputs_dir, output_dir = get_input_and_output_dirs(source_dir)
            os.makedirs(output_dir, exist_ok=True)

            log_fp = os.path.join(
                output_dir,
                f"{timestamp_str}_log_capture_and_plot_ladders.log")
            logger, curr_source_fhandler = set_up_logger(log_fp, config)

            per_sample_dfs_by_date = capture_qpcr_data(inputs_dir, config)
            concatted_per_sample_df = combine_all_datesets(
                per_sample_dfs_by_date)

            assay_dicts_list = config[INPUT_FILES_KEY][ASSAYS_KEY]
            for curr_assay_dict in assay_dicts_list:
                curr_ladders_fps, curr_plots_fps = \
                    generate_ladder_and_curve_plots_for_assay(
                        concatted_per_sample_df, curr_assay_dict,
                        timestamp_str, output_dir, config)
                ladder_fps.extend(curr_ladders_fps)
                plot_fps.extend(curr_plots_fps)
            # next assay

            curr_source_fhandler.close()
            logger.removeHandler(curr_source_fhandler)
            counter += 1
        # next source dir

        return ladder_fps, plot_fps
    except Exception:
        logging.exception("capture_and_save_ladders failed")
        raise

# TODO: maybe break this out to, e.g., output_formatted_data.py in case
#  we want to use this functionality again?
# def generate_model_training_files(
#         ladders_fps, runs_to_ignore, cq_conf_threshold, config):
#     RUN_COUNTER_KEY = "RunCounter_i"
#     RUNG_COUNTER_KEY = "CopyNumberCounter_j"
#     REP_COUNTER_KEY = "ReplicateCounter_k"
#     MEASUREMENT_KEY = "Cq_Yijk"
#     COPY_NUM_KEY = "CopyNumber_Xijk"
#
#     _, regression_bottom_step, _ = get_ladder_vals_from_config(config)
#
#     for curr_ladder_fp in ladders_fps:
#         curr_ladder_df = pandas.read_csv(curr_ladder_fp)
#
#         # TODO: note this hack won't work anymore--update if using again!
#         pieces = curr_ladder_fp.split("_")
#         joined_name_pieces = "_".join(pieces[-2:])
#         curr_target_name = joined_name_pieces.replace(".csv", "")
#
#         cq_conf_key = get_cq_conf_key_for_target(
#             curr_target_name, config)
#         run_mask = ~curr_ladder_df[RUN_NAME_KEY].isin(runs_to_ignore)
#         target_mask = curr_ladder_df[curr_target_name].notna()
#         conf_mask = curr_ladder_df[cq_conf_key] > cq_conf_threshold
#         reg_worthy_mask = \
#             curr_ladder_df[LADDER_STEP_KEY] >= regression_bottom_step
#         full_mask = run_mask & target_mask & conf_mask & reg_worthy_mask
#
#         working_df = curr_ladder_df.loc[full_mask, :].copy()
#         working_df.sort_values(
#             by=[RUN_NAME_KEY, LADDER_STEP_KEY, REPLICATE_KEY],
#             inplace=True)
#
#         # NB: ngroup is zero-based and we want 1-based
#         working_df[RUN_COUNTER_KEY] = \
#             working_df.groupby(RUN_NAME_KEY).ngroup() + 1
#         working_df[RUNG_COUNTER_KEY] = \
#             working_df.groupby(LADDER_STEP_KEY).ngroup() + 1
#         working_df.rename(columns={REPLICATE_KEY: REP_COUNTER_KEY,
#                                    curr_target_name: MEASUREMENT_KEY,
#                                    LADDER_STEP_KEY: COPY_NUM_KEY},
#                           inplace=True)
#
#         output_fp = curr_ladder_fp.replace(".csv", "_model_data.csv")
#         output_df = working_df.loc[:, [RUN_COUNTER_KEY,
#                                        RUNG_COUNTER_KEY,
#                                        REP_COUNTER_KEY,
#                                        MEASUREMENT_KEY,
#                                        COPY_NUM_KEY]].copy()
#         # output_df.sort_values(by=[RUN_COUNTER_KEY, RUNG_COUNTER_KEY,
#         #                           REP_COUNTER_KEY], inplace=True)
#         output_df.to_csv(output_fp, index=False)
#     # next ladder file


def generate_ladder_and_curve_plots_for_assay(
        cum_per_sample_df, assay_dict, timestamp_str, output_dir, config):
    ladders_output_fps = []
    assay_name = assay_dict[ASSAY_NAME_KEY]

    # get unique dates in the dataset
    sorted_unique_dataset_id_tuples = \
        get_unique_sorted_date_tuples(cum_per_sample_df)

    curve_pdf_fps_list = []
    curves_pdf_fp = f"{output_dir}/{timestamp_str}_"\
                    f"ladder_curves_{assay_name}.pdf"
    curves_pdf = None
    try:
        assay_targets = assay_dict[TARGETS_KEY]
        for curr_target_dict in assay_targets:
            target_name = curr_target_dict[TARGET_NAME_KEY]

            # get all the ladder values for that target across days
            curr_ladder_df = extract_ladders_for_target(
                cum_per_sample_df, curr_target_dict, config)

            if curr_ladder_df is None:
                logger = logging.getLogger()
                logger.info(f"No ladder data for {target_name}")
            else:
                if curves_pdf is None:
                    curves_pdf = PdfPages(curves_pdf_fp)
                    curve_pdf_fps_list.append(curves_pdf_fp)

                curr_ladders_output_fp = os.path.join(
                    output_dir, f"{timestamp_str}_ladder_{target_name}.csv")
                curr_ladder_df.to_csv(curr_ladders_output_fp, index=False)
                ladders_output_fps.append(curr_ladders_output_fp)

                curr_intercepts_per_run = None
                if target_name in config[STD_CURVE_KEY]:
                    curr_intercepts_per_run = \
                        get_mixed_model_intercepts_per_run_name(
                            curr_ladder_df, target_name, config)

                # generate the plot for this target
                _generate_std_curve_and_ladders_plot(
                    curr_ladder_df, target_name, sorted_unique_dataset_id_tuples,
                    curves_pdf, config, curr_intercepts_per_run)
            # endif ladders don't/do exist for this target
        # next target for assay
    finally:
        if curves_pdf is not None:
            curves_pdf.close()

    return ladders_output_fps, curve_pdf_fps_list


def get_single_curve_dicts_by_run(target_ladders_df, target_name, config,
                                  intercepts_per_run):
    single_curve_dict_per_run = {}
    mixed_curve_dict_per_run = {}

    sorted_unique_dataset_id_tuples = \
        get_unique_sorted_date_tuples(target_ladders_df)

    # get the target's standard curve
    target_curve_dict = config[STD_CURVE_KEY].get(target_name)
    if not target_curve_dict:
        logger = logging.getLogger()
        logger.info(f"No standard curve for {target_name}")
        target_curve_dict = {USE_LN_KEY: False}

    # for each date in the input dataset (NOT just the ladder set)
    for curr_dataset_id_tuple in sorted_unique_dataset_id_tuples:
        date_dataset_name = getattr(curr_dataset_id_tuple, RUN_NAME_KEY)
        ladder_for_date_df = _extract_target_ladder_for_date(
            target_ladders_df, curr_dataset_id_tuple, target_name)

        if ladder_for_date_df is not None:
            single_curve_dict, mixed_curve_dict = _get_curve_dicts_for_date(
                ladder_for_date_df, date_dataset_name, target_name,
                target_curve_dict, config, intercepts_per_run)

            single_curve_dict_per_run[date_dataset_name] = single_curve_dict
            mixed_curve_dict_per_run[date_dataset_name] = mixed_curve_dict
        # endif
    # next run/date

    return single_curve_dict_per_run, mixed_curve_dict_per_run


def get_mixed_model_intercepts_per_run_name(
        curr_ladder_df, curr_target_name, config):

    _, _, intercept_ladder_level = get_ladder_vals_from_config(config)

    curr_intercept_level_mask = \
        curr_ladder_df[LADDER_STEP_KEY] == intercept_ladder_level
    if not curr_intercept_level_mask.any():
        raise ValueError("No ladder values found at the intercept level")

    curr_target_cq_conf_key = get_cq_conf_key_for_target(
        curr_target_name, config)
    min_acceptable_cq_conf = \
        config[ANALYSIS_SETTINGS_KEY][MIN_ACCEPTABLE_CQ_CONF_KEY]
    cq_conf_passes_mask = \
        curr_ladder_df[curr_target_cq_conf_key] >= min_acceptable_cq_conf
    curr_combined_mask = curr_intercept_level_mask & cq_conf_passes_mask
    curr_intercept_ladder_df = \
        curr_ladder_df.loc[curr_combined_mask, :].copy()

    if len(curr_intercept_ladder_df) == 0:
        raise ValueError(f"Ladder measurements for intercept level failed Cq confidence")

    # TODO: should I put in some sort of check for there being
    #  at least some number for each run name?

    # NB: pivot_table, by default, averages the values for duplicate
    # indices, which is the desired behavior here.
    temp_df = curr_intercept_ladder_df.pivot_table(
        index=RUN_NAME_KEY, values=curr_target_name, dropna=False)

    target_curve_dict = config[STD_CURVE_KEY][curr_target_name]
    shared_slope = target_curve_dict[SLOPE_KEY]
    is_natural_log = target_curve_dict[USE_LN_KEY]
    if is_natural_log:
        exponent_term = numpy.log(intercept_ladder_level)
    else:  # is "common" log--ie., base 10
        exponent_term = numpy.log10(intercept_ladder_level)

    intercept_by_run_name = \
        temp_df[curr_target_name] - shared_slope * exponent_term
    return intercept_by_run_name


def extract_ladders_for_target(per_sample_df, target_dict, config):
    logger = logging.getLogger()

    target_name = target_dict[TARGET_NAME_KEY]
    if target_name not in per_sample_df.columns:
        logger.info(f"No data found for '{target_name}'")
        return None

    ladder_prefixes = target_dict.get(TARGET_LADDER_PREFIX_KEY)
    if not ladder_prefixes:
        logger.info(f"No ladder prefix(s) for target '{target_name}'")
        return None

    ladders_w_matching_prefixes_df = _extract_qpcr_ladder(
        per_sample_df, ladder_prefixes, config)
    if ladders_w_matching_prefixes_df is None:
        logger.info(f"No ladder values for {target_name}")
        return None

    # any records for which this target wasn't measured will be na;
    # exclude those!
    target_was_measured_mask = \
        ladders_w_matching_prefixes_df[target_name].notna()
    target_ladder_df = \
        ladders_w_matching_prefixes_df.loc[target_was_measured_mask, :].copy()

    # now that the "not-measured" na's are gone,
    # replace any "measured but not found" values with na,
    # so they are left out of calculations
    computable_target_ladder_df = convert_measured_but_not_found_to_nan(
        target_ladder_df)

    return computable_target_ladder_df


def _extract_qpcr_ladder(qpcr_sample_df, ladder_prefixes, config):
    # get all the entries from the qpcr data that start with the prefix(s);
    # if no ladder data found, return None
    ladder_pref_mask = qpcr_sample_df[INTERNAL_SAMPLE_KEY].str.startswith(
        tuple(ladder_prefixes))
    if not ladder_pref_mask.any():
        return None

    # get df of ladder records; add level for each. May be multiples per level
    ladder_vals_df = qpcr_sample_df.loc[ladder_pref_mask, :].copy()
    ladder_vals_df.sort_values(INTERNAL_SAMPLE_KEY, inplace=True)
    ladder_vals_df.reset_index(inplace=True, drop=True)

    ladder_info = ladder_vals_df.apply(
        lambda row: _extract_ladder_info_from_sample_name(
            row[INTERNAL_SAMPLE_KEY], ladder_prefixes, config),
        axis=1)

    # split column of lists into new columns then join back to orig df
    split = pandas.DataFrame(
        ladder_info.to_list(), columns=[LADDER_STEP_KEY, REPLICATE_KEY])
    ladder_vals_df = pandas.concat(
        [ladder_vals_df, split], axis=1)

    return ladder_vals_df


def _extract_ladder_info_from_sample_name(sample_name, ladder_prefixes, config):
    # first, remove the ladder prefix
    for curr_prefix in ladder_prefixes:
        sample_name = sample_name.replace(curr_prefix, "")

    # then split on replicate separator
    pieces = sample_name.split(REPLICATE_SEPARATOR)

    # take first piece and  convert that to a number
    ladder_level_str = pieces[0]
    # TODO: take out this temporary munge of past mis-labeled levels
    ladder_level_str = ladder_level_str.replace("6", "4")
    ladder_level = float(ladder_level_str)

    # then check if that number is in expected ladder levels and shout if not
    ladder_levels, _, _ = get_ladder_vals_from_config(config)
    if ladder_level not in ladder_levels:
        raise ValueError(f"Unexpected ladder level '{ladder_level}'")

    # take second piece (replicate num) and convert that to an integer
    ladder_replicate = 1  # default
    if len(pieces) > 1:
        ladder_replicate_str = pieces[1]
        ladder_replicate = int(ladder_replicate_str)

    return [ladder_level, ladder_replicate]


def get_ladder_vals_from_config(config):
    ladder_dict = config[INPUT_FILES_KEY][LADDERS_KEY]
    ladder_bottom = ladder_dict[LADDER_BOTTOM_KEY]
    ladder_steps = ladder_dict[LADDER_STEPS_KEY]
    ladder_vals = [ladder_bottom * (10 ** x) for x in range(0, ladder_steps)]

    regression_bottom_step = ladder_dict[REGRESSION_BOTTOM_STEP_KEY]
    regression_bottom_val = ladder_bottom * (10 ** regression_bottom_step)

    intercept_val = None
    intercept_step = ladder_dict.get(INTERCEPT_LADDER_STEP_KEY)
    if intercept_step is not None:
        intercept_val = ladder_bottom * (10 ** intercept_step)

    return ladder_vals, regression_bottom_val, intercept_val


def _generate_std_curve_and_ladders_plot(
        target_ladders_df, target_name, sorted_unique_dataset_id_tuples,
        curves_pdf, config, intercepts_per_run=None):

    # set up the graph object and add the target's standard curve
    target_curve_dict = config[STD_CURVE_KEY].get(target_name)
    if not target_curve_dict:
        logger = logging.getLogger()
        logger.info(f"No standard curve for {target_name}")
        # TODO: this seems pretty arbitrary; should these default values be
        #  specified in the config?
        target_curve_dict = {USE_LN_KEY: False,
                             REACTION_VOL_PER_L_KEY: 1000000}
    curve_plot = _set_up_std_curve_plot(
        target_curve_dict, target_name, config)

    # for each date in the input dataset (NOT just the ladder set)
    for curr_dataset_id_tuple in sorted_unique_dataset_id_tuples:
        date_dataset_name = getattr(curr_dataset_id_tuple, RUN_NAME_KEY)
        ladder_for_date_df = _extract_target_ladder_for_date(
            target_ladders_df, curr_dataset_id_tuple, target_name)

        if ladder_for_date_df is not None:
            single_curve_dict, mixed_curve_dict = _get_curve_dicts_for_date(
                ladder_for_date_df, date_dataset_name, target_name,
                target_curve_dict, config, intercepts_per_run)

            # TODO: if we decide to keep using mixed model curves, we will
            #  probably want to figure out a more elegant way than adding a
            #  temporary config key
            if not config[INCLUDE_MIXED_KEY]:
                mixed_curve_dict = None

            curve_plot = _add_date_dataset_to_plot(
                curve_plot, ladder_for_date_df, date_dataset_name, target_name,
                single_curve_dict, mixed_curve_dict)

    _save_curve_plot(curve_plot, curves_pdf)


def _set_up_std_curve_plot(target_curve_dict, title, config):
    # The pylab figure manager will be bypassed in this instance.
    # This means that `fig` will be garbage collected as you'd expect.
    fig = Figure()
    ax = fig.subplots()

    ladder_copies_per_reaction, fit_bottom_val, _ = \
        get_ladder_vals_from_config(config)

    ax, base = _plot_regression_line(
        target_curve_dict, ladder_copies_per_reaction, ax)

    # Label plot
    ax.set_xscale("log", base=base)
    ax.set_xlabel('gene copies/uL')
    ax.set_ylabel('Cq')
    ax.set_title(f"{title} ladders and curves (fit down to "
                 f"{fit_bottom_val} copies/uL)")

    ax.yaxis.set_major_locator(MultipleLocator(3))
    ax.yaxis.set_minor_locator(MultipleLocator(1))
    ax.grid(visible=True, which="both", axis='y')

    return ax


def _plot_regression_line(target_curve_dict, ladder_copies_per_reaction, ax):
    # Plot ladder regression line
    if target_curve_dict[USE_LN_KEY]:
        xseq = numpy.log(ladder_copies_per_reaction)
        base = numpy.e
        label_txt = "ln"
    else:
        xseq = numpy.log10(ladder_copies_per_reaction)
        base = 10
        label_txt = "log10"

    if SLOPE_KEY in target_curve_dict:
        label = f"std curve ({target_curve_dict[SLOPE_KEY]} {label_txt}(x) + " \
                f"{target_curve_dict[INTERCEPT_KEY]}) "
        ax.plot(ladder_copies_per_reaction,
                target_curve_dict[INTERCEPT_KEY] + target_curve_dict[SLOPE_KEY] * xseq,
                lw=1.5, label=label, color="black")

    return ax, base


def _add_date_dataset_to_plot(
        working_plot, ladder_for_date_df, date_dataset_name,
        target_name, single_curve_dict, mixed_curve_dict):

    all_ladder_copies_per_reaction, all_ladder_cqs = \
        _extract_ladder_copies_and_cq_lists(
            ladder_for_date_df, target_name)

    working_plot = _add_ladder_points_and_curve_to_plot(
        working_plot, single_curve_dict, all_ladder_copies_per_reaction,
        all_ladder_cqs, date_dataset_name, mixed_curve_dict)

    return working_plot


def _get_curve_dicts_for_date(ladder_for_date_df, date_dataset_name,
                              target_name, target_curve_dict, config,
                              intercepts_per_run=None):

    if ladder_for_date_df is not None:
        _, fit_bottom_val, _ = get_ladder_vals_from_config(config)
        has_enough_ladder_copies_for_fit_mask = \
            ladder_for_date_df[LADDER_STEP_KEY] >= fit_bottom_val
        fit_ladder_copies_per_reaction, fit_ladder_cqs = \
            _extract_ladder_copies_and_cq_lists(
                ladder_for_date_df, target_name,
                has_enough_ladder_copies_for_fit_mask)

        single_model_curve_dict = _calc_curve_from_ladder(
            fit_ladder_copies_per_reaction, fit_ladder_cqs,
            target_curve_dict[USE_LN_KEY])
        single_model_curve_dict[REACTION_VOL_PER_L_KEY] = \
            target_curve_dict[REACTION_VOL_PER_L_KEY]
        single_model_curve_dict[CURVE_NAME_KEY] = \
            f"{SINGLE_CURVE_VAL}_{CURVE_MODEL_KEY}"

        mixed_model_curve_dict = None
        if intercepts_per_run is not None:
            mixed_model_curve_dict = target_curve_dict.copy()
            mixed_model_curve_dict[INTERCEPT_KEY] = \
                intercepts_per_run[date_dataset_name]
            mixed_model_curve_dict[CURVE_NAME_KEY] = \
                f"{MIXED_CURVE_VAL}_{CURVE_MODEL_KEY}"

        return single_model_curve_dict, mixed_model_curve_dict


def _extract_target_ladder_for_date(
        target_ladder_df, dataset_id_tuple, target_name):

    dataset_date = getattr(dataset_id_tuple, RUN_DATE_KEY)
    dataset_date_mask = \
        target_ladder_df[RUN_DATE_KEY] == dataset_date
    dataset_date_df = target_ladder_df.loc[dataset_date_mask, :]
    if len(dataset_date_df) == 0 or dataset_date_df[target_name].isna().all():
        dataset_name = getattr(dataset_id_tuple, RUN_NAME_KEY)
        logger = logging.getLogger()
        logger.info(f"No ladder values for {target_name} in {dataset_name}")
        return None

    # TODO: should I put in check to make sure it has ALL the expected levels?
    return dataset_date_df


def _extract_ladder_copies_and_cq_lists(ladder_df, target_name, mask=None):
    working_df = ladder_df

    if mask is not None:
        working_df = ladder_df.loc[mask, :]

    ladder_copies_per_reaction = working_df[LADDER_STEP_KEY].tolist()
    ladder_raw_cqs = working_df[target_name].tolist()

    return ladder_copies_per_reaction, ladder_raw_cqs


def _calc_curve_from_ladder(ladder_copies_per_reaction, ladder_raw_cqs,
                            use_natural_log=True):
    # calculate coefficients of the regression line (1st degree polynomial)
    # for the ladder--lns or log10s of ladder copies/reaction are the
    # independent variables, ladder Cqs are the dependent ones
    result = {}

    raw_cqs_array = numpy.array(ladder_raw_cqs)
    copies_per_reaction_array = numpy.array(ladder_copies_per_reaction)
    # NB: np.log is ln, whereas np.log10 is standard base 10 log
    if use_natural_log:
        log_cp_per_rxn_array = numpy.log(copies_per_reaction_array)
    else:
        log_cp_per_rxn_array = numpy.log10(copies_per_reaction_array)

    # include only indices that are *not* NaN in either array
    non_nan_idx = numpy.isfinite(log_cp_per_rxn_array) & numpy.isfinite(
        raw_cqs_array)
    non_nan_log_cp_per_rxn = log_cp_per_rxn_array[non_nan_idx]
    non_nan_raw_cqs = raw_cqs_array[non_nan_idx]
    coeffs = numpy.polyfit(non_nan_log_cp_per_rxn, non_nan_raw_cqs, deg=1)

    # r
    correlation = numpy.corrcoef(non_nan_log_cp_per_rxn, non_nan_raw_cqs)[0, 1]
    # r-squared
    determination = correlation**2

    result[SLOPE_KEY] = coeffs[0]
    result[INTERCEPT_KEY] = coeffs[1]
    result[R_SQUARED_KEY] = determination
    result[USE_LN_KEY] = use_natural_log
    return result


def _add_ladder_points_and_curve_to_plot(
        ax, ladder_curve_dict, ladder_copies_per_reaction,
        ladder_raw_cqs, ladder_title, mixed_model_curve_dict):

    color = next(ax._get_lines.prop_cycler)['color']

    # Make scatter plot of ladder data points
    ax.scatter(ladder_copies_per_reaction, ladder_raw_cqs,
               s=8, alpha=0.7, label=f"{ladder_title} ladder", color=color)

    ax = _add_curve_to_plot(ax, ladder_curve_dict, ladder_copies_per_reaction, color=color)

    if mixed_model_curve_dict is not None:
        ax = _add_curve_to_plot(ax, mixed_model_curve_dict,
                                ladder_copies_per_reaction, use_dotted=True, color=color)

    return ax


def _add_curve_to_plot(ax, ladder_curve_dict, ladder_copies_per_reaction,
                       use_dotted=False, color=None):
    # Plot ladder regression line
    if ladder_curve_dict[USE_LN_KEY]:
        xseq = numpy.log(ladder_copies_per_reaction)
        label_txt = "ln"
    else:
        xseq = numpy.log10(ladder_copies_per_reaction)
        label_txt = "log10"

    linestyle = "solid"
    if use_dotted:
        linestyle = "dashed"

    r_squared_text = ""
    if R_SQUARED_KEY in ladder_curve_dict:
        r_squared_text = f", R^2 = {round(ladder_curve_dict[R_SQUARED_KEY], 4)}"

    # NB: took out ladder title from line label (because they are long and
    # messy) so now line is associated with scatter just by shared color.
    label = f"{round(ladder_curve_dict[SLOPE_KEY], 3)} {label_txt}(x) " \
            f"+ {round(ladder_curve_dict[INTERCEPT_KEY], 3)}{r_squared_text}"
    ax.plot(ladder_copies_per_reaction,
            ladder_curve_dict[INTERCEPT_KEY] + ladder_curve_dict[SLOPE_KEY] * xseq,
            lw=1, label=label, linestyle=linestyle, color=color)

    return ax


def _save_curve_plot(target_curve_plot, pdf_pages):
    curr_lgd = target_curve_plot.legend(loc=9, bbox_to_anchor=(0.5, -0.1))

    # The below is to ensure that figures are garbage-collected; see
    # https://stackoverflow.com/q/21884271 ,
    # https://stackoverflow.com/a/16337909 ,
    # https://stackoverflow.com/a/16337909#comment23407346_16337909
    # Note the need to create canvas even though not being used:
    # "It's a wart in the API; ideally the canvas would be initiated along
    # with the figure" (but it isn't)
    FigureCanvas(target_curve_plot.figure)
    pdf_pages.savefig(target_curve_plot.figure, bbox_extra_artists=(curr_lgd,),
                      bbox_inches='tight')
