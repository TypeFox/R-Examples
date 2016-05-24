
trace_length_case <- function(eventlog) {
	stop_eventlog(eventlog)

	case_classifier <- case_id(eventlog)
	colnames(eventlog)[colnames(eventlog) == activity_instance_id(eventlog)] <- "activity_instance_classifier"

	r<- eventlog %>% group_by_(case_classifier) %>% summarize(trace_length = n_distinct(activity_instance_classifier))

	return(r)
}
