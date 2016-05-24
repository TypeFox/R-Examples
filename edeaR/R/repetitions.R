
#' @title Metric:  Repetitions
#'
#' @description  Provides summuary statistics on the number of repetitions, at the level of activity types, traces, cases and the eventlog.
#'
#' @param eventlog The event log to be used. An object of class
#' \code{eventlog}.
#'
#' @param level_of_analysis At which level the analysis of repetitions should be performed: log, case, trace or activity.
#'
#'
#' @export repetitions


repetitions <- function(eventlog,
						level_of_analysis){

	stop_eventlog(eventlog)

	if(!(level_of_analysis %in% c("log","trace","case","activity")))
		stop("Level of analysis should be one of the following: log, trace, case, activity.")

	if(level_of_analysis == "log"){
		return(repetitions_log(eventlog = eventlog))
	}
	else if(level_of_analysis == "trace") {
		return(repetitions_trace(eventlog = eventlog))
	}
	else if(level_of_analysis == "case")
		return(repetitions_case(eventlog = eventlog))
	else
		return(repetitions_activity(eventlog = eventlog))
}
