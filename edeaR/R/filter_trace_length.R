#' @title Filter: Trace length percentile
#'
#' @description Filters cases on length, using a percentile threshold.
#'
#' @param eventlog The event log to be used. An object of class
#' \code{eventlog}.
#'
#' @param lower_threshold The lower trace length threshold.
#' When \code{reverse} is FALSE, all traces with a lower frequency are discarded.
#'
#' @param upper_threshold The upper trace length threshold.
#' When \code{reverse} is FALSE, all traces with a lo frequency are discarded.
#'
#' @param percentile_cut_off Alternatively to providing thresholds, a percentile cut off can be provided.
#' A percentile cut off value of 0.9 will return the 90\% shortest cases.
#' When \code{reverse} is set to TRUE, it will return the 10\% longest cases.
#'
#' @param reverse A logical parameter depicting whether the selection should be reversed.
#'
#'
#' @export filter_trace_length
#'


filter_trace_length <- function(eventlog,
								lower_threshold = NULL,
								upper_threshold = NULL,
								percentile_cut_off = NULL,
								reverse = F) {

	stop_eventlog(eventlog)

	if(is.null(lower_threshold) & is.null(upper_threshold) & is.null(percentile_cut_off))
		stop("At least one threshold or a percentile cut off must be provided.")


	if((!is.null(lower_threshold) & !is.null(percentile_cut_off)) | (!is.null(upper_threshold) & !is.null(percentile_cut_off)))
		stop("Cannot filter on both thresholds and percentile cut off simultaneously.")


	if(!is.null(percentile_cut_off))
		return(filter_trace_length_percentile(eventlog,
											  percentile_cut_off = percentile_cut_off,
											  reverse = reverse ))
	else
		return(filter_trace_length_threshold(eventlog,
											 lower_threshold = lower_threshold,
											 upper_threshold = upper_threshold,
											 reverse = reverse))
}
