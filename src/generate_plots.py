import pandas
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.dates as mdates
from matplotlib.ticker import MultipleLocator, StrMethodFormatter
from src.shared_symbols import *
from src.utils import extract_target_and_control_names, get_assay_dict_by_name
from src.capture_qpcr_data import get_unique_sorted_date_tuples

FONT_SIZE = 4
DATE_TITLE = "Collection Date"
TABLEAU_COLORS_MOD = ["red", "blue", "orange", "green", "purple", "brown",
                      "pink", "gray", "olive", "cyan"]


def generate_plots(processing_run_name, dirnames, sample_dfs, base_dfs,
                   output_dir, config,
                   earliest_date=None, do_zoomed=False):
    output_assay_dict = config[OUTPUT_FILES_KEY][OUTPUT_ASSAY_KEY]
    assay_name = output_assay_dict[ASSAY_NAME_KEY]
    assay_dict = get_assay_dict_by_name(assay_name, config)

    cqs_pdf = PdfPages(f"{output_dir}/{processing_run_name}_cqs.pdf")
    concs_pdf = PdfPages(f"{output_dir}/{processing_run_name}_concs.pdf")
    zoomed_cqs_pdf = None
    zoomed_concs_pdf = None
    if do_zoomed:
        zoomed_cqs_pdf = PdfPages(f"{output_dir}/{processing_run_name}_zoomed_cqs.pdf")
        zoomed_concs_pdf = PdfPages(f"{output_dir}/{processing_run_name}_zoomed_concs.pdf")

    try:
        for curr_idx in range(0, len(dirnames)):
            # TODO: figure out how to incorporate dirname
            curr_dirname = dirnames[curr_idx]
            curr_sample_df = sample_dfs[curr_idx]
            curr_base_df = base_dfs[curr_idx]

            if curr_sample_df is not None and len(curr_sample_df) > 0:
                make_and_save_qc_plots(
                    curr_sample_df, curr_base_df, assay_dict, cqs_pdf, concs_pdf,
                    zoomed_cqs_pdf, zoomed_concs_pdf, earliest_date)
        # next dirname
    finally:
        cqs_pdf.close()
        concs_pdf.close()
        if do_zoomed:
            zoomed_cqs_pdf.close()
            zoomed_concs_pdf.close()


def make_and_save_qc_plots(per_sample_df, per_base_df, assay_dict, cqs_pdf,
                           concs_pdf, zoomed_cqs_pdf, zoomed_concs_pdf,
                           earliest_date):
    unique_loc_prefixes, unique_run_names_ordered_by_desc_run_date = \
        _get_unique_locs_and_run_names(per_sample_df)

    for curr_loc_prefix in unique_loc_prefixes:
        curr_per_base_loc_subset_df = _get_location_subset_df(
            per_base_df, curr_loc_prefix, earliest_date)

        curr_cq_title = f"Raw Cqs for {curr_loc_prefix}"
        curr_conc_title = f"Viral Load for {curr_loc_prefix}"
        curr_cq_plot = None
        curr_zoomed_cq_plot = None
        curr_norm_conc_plot = None
        curr_zoomed_norm_conc_plot = None

        curr_color_index = 0
        for curr_run_name in unique_run_names_ordered_by_desc_run_date:
            if curr_color_index == len(TABLEAU_COLORS_MOD):
                curr_color_index = curr_color_index - len(TABLEAU_COLORS_MOD)
            run_color = TABLEAU_COLORS_MOD[curr_color_index]

            curr_per_sample_subset_df = _get_run_and_location_subset_df(
                per_sample_df, curr_run_name, curr_loc_prefix, earliest_date)
            if len(curr_per_sample_subset_df) > 0:
                curr_cq_plot = _make_or_add_to_plot(
                    curr_cq_plot, curr_per_sample_subset_df, curr_run_name,
                    run_color, curr_cq_title, DATE_TITLE, "Cq", assay_dict)

                if zoomed_cqs_pdf:
                    curr_zoomed_cq_plot = _make_or_add_to_plot(
                        curr_zoomed_cq_plot, curr_per_sample_subset_df, curr_run_name,
                        run_color, curr_cq_title, DATE_TITLE, "Cq", assay_dict,
                        zoomed=True)

                curr_norm_conc_plot = _make_or_add_to_plot(
                    curr_norm_conc_plot, curr_per_sample_subset_df,
                    curr_run_name, run_color, curr_conc_title, DATE_TITLE,
                    "Gene Copies per Liter")

                if zoomed_concs_pdf:
                    curr_zoomed_norm_conc_plot = _make_or_add_to_plot(
                        curr_zoomed_norm_conc_plot, curr_per_sample_subset_df,
                        curr_run_name, run_color, curr_conc_title, DATE_TITLE,
                        "Gene Copies per Liter", zoomed=True)
                # endif no plot for this loc yet

                if len(curr_per_base_loc_subset_df) == 0:
                    raise ValueError(
                        f"No per-base entries for location "
                        f"'{curr_loc_prefix}' that has per-sample entries")
                # endif there aren't any per-base records
            # end if there are records for this run/location

            curr_color_index += 1
        # next run

        curr_norm_conc_plot = _add_concs_scatter(
            curr_per_base_loc_subset_df, "Averaged",
            curr_norm_conc_plot, "black",
            fails_qc_key=IGNORE_KEY,
            make_big=True)

        if zoomed_concs_pdf:
            curr_zoomed_norm_conc_plot = _add_concs_scatter(
                curr_per_base_loc_subset_df, "Averaged",
                curr_zoomed_norm_conc_plot, "black",
                fails_qc_key=IGNORE_KEY,
                make_big=True)

        _save_plot_to_multipdf(cqs_pdf, curr_cq_plot)
        if zoomed_cqs_pdf:
            _save_plot_to_multipdf(zoomed_cqs_pdf, curr_zoomed_cq_plot)
        _save_plot_to_multipdf(concs_pdf, curr_norm_conc_plot)
        if zoomed_concs_pdf:
            _save_plot_to_multipdf(zoomed_concs_pdf, curr_zoomed_norm_conc_plot)
    # next loc prefix


def _get_unique_locs_and_run_names(per_sample_df):
    unique_loc_prefixes = pandas.unique(per_sample_df.loc[:, LOCATION_KEY])

    sorted_unique_dataset_id_tuples = get_unique_sorted_date_tuples(
        per_sample_df, sort_asc=False)
    unique_run_names_ordered_by_desc_run_date = \
        [getattr(x, RUN_NAME_KEY) for x in sorted_unique_dataset_id_tuples]

    return unique_loc_prefixes, unique_run_names_ordered_by_desc_run_date


def _make_or_add_to_plot(a_plot, a_df, run_name, color,
                         title, x_title, y_title, assay_dict=None,
                         zoomed=False):
    x_units = None
    y_units = None
    y_min = None
    y_max = None
    if a_plot is None:
        if zoomed:
            x_units = 1
            if not assay_dict:
                y_max = 10000000
                y_min = -1
                y_units = 1000000
            else:
                y_min = 19
                y_max = 40
                y_units = 2

        a_plot = _set_up_plot(title, x_title, y_title,
                              x_units=x_units, y_units=y_units,
                              y_min=y_min, y_max=y_max)
    # endif no plot yet

    truncated_run_name = run_name.rsplit("_Results")[0]
    if assay_dict:
        a_plot = _add_cqs_scatter(
            a_df, truncated_run_name, a_plot, color, assay_dict)
    else:
        a_plot = _add_concs_scatter(a_df, truncated_run_name, a_plot, color)

    return a_plot


def _set_up_plot(title, x_axis_title, y_axis_title, x_units, y_units,
                 y_min, y_max):
    # The pylab figure manager will be bypassed in this instance.
    # This means that `fig` will be garbage collected as you'd expect.
    fig = Figure()
    ax = fig.subplots()

    # Label plot
    ax.set_xlabel(x_axis_title, fontsize=FONT_SIZE)
    ax.set_ylabel(y_axis_title, fontsize=FONT_SIZE)
    ax.set_title(title)

    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%Y'))
    ax.yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    if x_units:
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=x_units))
    if y_units:
        ax.yaxis.set_major_locator(MultipleLocator(y_units))
        ax.yaxis.set_minor_locator(MultipleLocator(y_units/2))
    if y_min and y_max:
        ax.set_ylim(y_min, y_max)

    ax.tick_params(axis='both', labelsize=FONT_SIZE)
#    ax.ticklabel_format(axis="y", style="plain")
    ax.grid(visible=True, which="both", axis='y')
    fig.autofmt_xdate(ha='left', rotation=-45)

    return ax


def _get_run_and_location_subset_df(input_df, run_name, location_prefix,
                                    earliest_date=None):
    loc_subset_df = _get_location_subset_df(
        input_df, location_prefix, earliest_date)
    subset_df = _get_run_subset_df(loc_subset_df, run_name)
    return subset_df


def _get_location_subset_df(input_df, location_prefix, earliest_date=None):
    loc_mask = input_df[LOCATION_KEY] == location_prefix
    if earliest_date is not None:
        date_mask = input_df[INTERNAL_COLLECTION_DATE_KEY] >= earliest_date
        loc_mask = loc_mask & date_mask
    subset_df = input_df.loc[loc_mask, :].copy()
    return subset_df


def _get_run_subset_df(input_df, run_name):
    run_mask = input_df[RUN_NAME_KEY] == run_name
    subset_df = input_df.loc[run_mask, :].copy()
    return subset_df


def _add_cqs_scatter(run_and_location_per_sample_df, run_name,
                     curr_graph, curr_color, assay_dict):

    target_name, control_name = extract_target_and_control_names(assay_dict)
    for curr_target_name in [target_name, control_name]:
        for curr_qc_fail_status in [False, True]:
            not_target = curr_target_name != target_name
            curr_graph = _add_per_date_scatter_w_qc_constraint(
                run_and_location_per_sample_df, curr_target_name,
                FAILS_SAMPLE_QC_KEY, curr_qc_fail_status, run_name,
                curr_color, None, not_target, curr_graph)
        # next qc fail status
    # next target

    return curr_graph


def _add_per_date_scatter_w_qc_constraint(
        input_df, vals_key, fails_qc_key, fails_qc_val, series_name, color,
        marker_size, use_alternate_markers, a_plot):
    qc_mask = input_df[fails_qc_key] == fails_qc_val
    has_date_mask = ~input_df[INTERNAL_COLLECTION_DATE_KEY].isna()
    combined_mask = qc_mask & has_date_mask
    temp_df = input_df.loc[combined_mask, :].copy()

    if len(temp_df) > 0:
        datetimes = (temp_df.loc[:, INTERNAL_COLLECTION_DATE_KEY]).tolist()
        dates = [d.date() for d in datetimes]
        vals = (temp_df.loc[:, vals_key]).tolist()
    else:
        # seems weird to plot an empty series, but it is necessary in order to
        # get the appropriate legend entry for this run.  Note that we only get
        # into this function if there is SOME data for the run, so here we are
        # missing one kind of data (e.g., passing PMMoV cqs) but we know that
        # there is other relevant data for this run that will be plotted (and
        # thus this run will need a label).
        dates = []
        vals = []

    series_title = f"{series_name} {vals_key}"
    a_plot = _add_scatter(a_plot, dates, vals, series_title, color,
                          marker_size,
                          use_ignore_markers=fails_qc_val,
                          use_alternate_markers=use_alternate_markers)
    return a_plot


def _add_scatter(ax, x_vals, y_vals, series_title, color, marker_size,
                 use_ignore_markers=False, use_alternate_markers=False):
    label = None
    marker_size = 4 if marker_size is None else marker_size
    if use_ignore_markers:
        marker_type = '_' if use_alternate_markers else "x"
    else:
        if use_alternate_markers:
            marker_type = "^"
        else:
            marker_type = "."
            label = series_title
        # endif this is the "real" data
    # endif use/don't use ignore markers

    # Make scatter plot of ladder data points
    ax.scatter(x_vals, y_vals, s=marker_size, linewidth=1, label=label,
               facecolor=color, edgecolor=None, marker=marker_type)
    return ax


def _add_concs_scatter(input_df, series_name, curr_graph, curr_color,
                       fails_qc_key=IGNORE_KEY, make_big=False):
    marker_size = 12 if make_big else None
    for curr_qc_fail_status in [False, True]:
        curr_graph = _add_per_date_scatter_w_qc_constraint(
            input_df, INTERNAL_NORM_CONC_KEY,
            fails_qc_key, curr_qc_fail_status, series_name,
            curr_color, marker_size, False, curr_graph)
    # next qc fail status

    return curr_graph


def _save_plot_to_multipdf(pdf_pages, a_plot):
    if a_plot is None:
        return

    # TODO if keeping this, come up w better naming/storage scheme
    curr_lgd = a_plot.legend(loc=9, bbox_to_anchor=(0.5, -0.2),
                             fontsize=FONT_SIZE)

    # The below is to ensure that figures are garbage-collected; see
    # https://stackoverflow.com/q/21884271 ,
    # https://stackoverflow.com/a/16337909 ,
    # https://stackoverflow.com/a/16337909#comment23407346_16337909
    # Note the need to create canvas even though not being used:
    # "It's a wart in the API; ideally the canvas would be initiated along
    # with the figure" (but it isn't)
    canvas = FigureCanvas(a_plot.figure)
    pdf_pages.savefig(a_plot.figure, bbox_extra_artists=(curr_lgd,),
                      bbox_inches='tight')
