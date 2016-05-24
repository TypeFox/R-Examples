#' @title Filter: Activity frequency
#'
#' @description Filters the log based on its most frequent activities, until a specific percentile cut off.
#'
#' @param eventlog The event log to be used. An object of class
#' \code{eventlog}.
#'
#' @param percentile_cut_off The target coverage of events
#' A percentile of 0.9 will return the most common activity types of the eventlog, which account for 90\% of the events.
#'
#' @param reverse A logical parameter depicting whether the selection should be reversed.
#'
#' @export filter_activity_frequency
#'
filter_activity_frequency <- function(eventlog,
									 percentile_cut_off = 0.8,
									 reverse = F){
	stop_eventlog(eventlog)

	act_freq <- mutate(activities(eventlog), r = percent_rank(-absolute_frequency))

	if(reverse == F)
		event_selection <- act_freq %>% filter(r <= percentile_cut_off)

	else
		event_selection <- act_freq %>% filter(r > percentile_cut_off)


	colnames(event_selection)[colnames(event_selection) == activity_id(eventlog)] <- "event_classifier"
	colnames(eventlog)[colnames(eventlog) == activity_id(eventlog)] <- "event_classifier"

	event_selection <- select(event_selection, event_classifier)

	output <- filter(eventlog, event_classifier %in% event_selection$event_classifier)

	colnames(output)[colnames(output)=="event_classifier"] <- activity_id(eventlog)

	output <- eventlog(output,
					   activity_id = activity_id(eventlog),
					   case_id = case_id(eventlog),
					   timestamp =timestamp(eventlog),
					   lifecycle_id = lifecycle_id(eventlog),
					   activity_instance_id = activity_instance_id(eventlog))

	return(output)
}
