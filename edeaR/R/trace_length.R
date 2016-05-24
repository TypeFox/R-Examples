#' @title Metric: Trace length
#'
#' @description Computes the length of each trace, in terms of the number of events, at the level of the eventlog or the level of a trace.
#' The relative numbers at trace level measure trace length compared to the average trace length of the top 80% cases, approximately.
#'
#'
#' @param eventlog The event log to be used. An object of class
#' \code{eventlog}.
#' @param level_of_analysis At which level the analysis of  trace_length should be performed: log, case or trace.

#' @export trace_length

trace_length <- function(eventlog,
						 level_of_analysis) {

	stop_eventlog(eventlog)


	if(!(level_of_analysis %in% c("log","case","trace")))
		stop("Level of analysis should be one of the following: log, case, trace.")

	if(level_of_analysis == "trace")
		return(trace_length_trace(eventlog = eventlog))
	else if (level_of_analysis == "case")
		return(trace_length_case(eventlog = eventlog))
	else
		return(trace_length_log(eventlog = eventlog))

}
