#' @title Metric: Size of selfloops
#'
#' @description Provides summary statistics on the sizes of selfloops at the level of activity types, cases, traces or log. A selfloop of size x refers to the occurrence of x consecutive events
#' of that activity type.
#'
#' @param eventlog The event log to be used. An object of class
#' \code{eventlog}.
#'
#' @param level_of_analysis At which level the analysis of selfloops should be performed: log, case, trace or activity.
#'
#'
#' @export size_of_selfloops

size_of_selfloops <- function(eventlog,
							  level_of_analysis){

	stop_eventlog(eventlog)


	if(!(level_of_analysis %in% c("trace","activity","case", "log")))
		stop("Level of analysis should be one of the following: log, trace, activity.")


	if(level_of_analysis == "trace") {
		return(size_of_selfloops_trace(eventlog = eventlog))
	}
	else if(level_of_analysis == "activity"){
		return(size_of_selfloops_activity(eventlog = eventlog))
	}
	else if(level_of_analysis == "case") {
		return(size_of_selfloops_case(eventlog = eventlog))
	}
	else {
		return(size_of_selfloops_log(eventlog = eventlog))
	}
}
