

processing_time_case <- function(eventlog,
								 units = "days") {
	stop_eventlog(eventlog)

	colnames(eventlog)[colnames(eventlog) == timestamp(eventlog)] <- "timestamp_classifier"

	case_classifier <- case_id(eventlog)
	event_classifier <- activity_id(eventlog)
	activity_instance_classifier <- activity_instance_id(eventlog)

	e <- eventlog %>%
		group_by_(case_classifier, activity_instance_classifier) %>%
		summarize(s = min(timestamp_classifier), e = max(timestamp_classifier)) %>%
		mutate(dur = as.double(e - s , units = units)) %>%
		summarize(tot_dur = sum(dur)) %>%
		select_(case_classifier, "tot_dur") %>%
		rename(processing_time = tot_dur) %>%
		arrange(desc(processing_time))

	return(e)

}
