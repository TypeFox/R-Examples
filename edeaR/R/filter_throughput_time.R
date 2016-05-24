#' @title Filter: Throughput Time
#'
#' @description Filters cases based on their throughput time.
#'
#' @param eventlog The event log to be used. An object of class
#' \code{eventlog}.
#'
#' @param lower_threshold The lower duration threshold, specified in number of days.
#' When \code{reverse} is FALSE, all cases with a lower duration are discarded.
#'
#' @param upper_threshold The upper duration threshold, specified in number of days.
#' When \code{reverse} is FALSE, all cases with a higher duration are discarded.
#'
#' @param percentile_cut_off Alternatively to providing thresholds, a percentile cut off can be provided.
#' A percentile cut off value of 0.9 will return the 90\% shortest cases.
#' When \code{reverse} is set to TRUE, it will return the 10\% longest cases.
#'
#' @param reverse A logical parameter depicting whether the selection should be reversed.
#'
#'
#' @export filter_throughput_time
#'

filter_throughput_time <- function(eventlog,
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
		return(filter_throughput_time_percentile(eventlog,
												 percentile_cut_off = percentile_cut_off,
												 reverse = reverse ))
	else
		return(filter_throughput_time_threshold(eventlog,
												lower_threshold = lower_threshold,
												upper_threshold = upper_threshold,
												reverse = reverse))
}
