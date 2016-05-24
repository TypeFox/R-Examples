end_activities_case <- function(eventlog) {

	stop_eventlog(eventlog)

	case_classifier <- case_id(eventlog)
	timestamp_classifier <- timestamp(eventlog)
	event_classifier <- activity_id(eventlog)
	colnames(eventlog)[colnames(eventlog) == activity_id(eventlog)] <- "event_classifier"

	r <- eventlog %>% group_by_(case_classifier) %>% arrange_(timestamp_classifier) %>% summarize(end_activity = last(event_classifier))

	return(r)

}
