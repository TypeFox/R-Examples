

throughput_time_case <- function(eventlog, units = "days") {
	stop_eventlog(eventlog)

	case_classifier <- case_id(eventlog)
	colnames(eventlog)[colnames(eventlog) == timestamp(eventlog)] <- "timestamp_classifier"

	e <- eventlog %>%
		group_by_(case_classifier) %>%
		summarize(s = min(timestamp_classifier), e = max(timestamp_classifier)) %>%
		mutate(throughput_time = as.double( e - s, units = units)) %>%
		arrange(throughput_time) %>% select(-s, -e)

	return(e)


}
