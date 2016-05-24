#' @title Metric: Processing time
#'
#' @description Provides summary statistics about the processing time of events on the level of activities, traces, cases or log.
#'
#' @param eventlog The event log to be used. An object of class
#' \code{eventlog}.
#'
#' @param units The time unit in which the throughput times should be reported.
#'
#' @param level_of_analysis At which level the analysis of processing times should be performed: log, trace, case or activity.
#'
#' @export processing_time

processing_time <- function(eventlog,
							level_of_analysis,
							units = "days"){

	stop_eventlog(eventlog)

	if(!(level_of_analysis %in% c("log","trace","case","activity")))
		stop("Level of analysis should be one of the following: log, trace, case,  activity.")

	if(level_of_analysis == "trace"){
		return(processing_time_trace(eventlog = eventlog, units = units))
	}
	else if(level_of_analysis == "log") {
		return(processing_time_log(eventlog = eventlog, units = units))
	}
	else if(level_of_analysis == "case") {
		return(processing_time_case(eventlog = eventlog, units = units))
	}
	else {
		return(processing_time_activity(eventlog = eventlog, units = units))
	}
}
