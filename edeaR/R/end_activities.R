#'@title Metric: End activities
#'
#'@description At log level, computes how many activity types occur as the last event in a case, both absolute and relative.
#'At activity level, shows the activities which occur as last, and how often.
#'The last event in a case is the one which completed the last.
#'
#' @param level_of_analysis At which level the analysis of end activities should be performed: log, case or activity.
#'
#' @param eventlog The event log to be used. An object of class
#' \code{eventlog}.
#'

#'@export end_activities


end_activities <- function(eventlog,
						   level_of_analysis) {
	stop_eventlog(eventlog)

	if(!(level_of_analysis %in% c("log","case", "activity")))
		stop("Level of analysis should be one of the following: log, case, activity.")

	if (level_of_analysis == "log")
		return(end_activities_log(eventlog = eventlog))
	else if (level_of_analysis == "case")
		return(end_activities_case(eventlog = eventlog))
	else
		return(end_activities_activity(eventlog = eventlog))
}
