size_of_selfloops_case <- function(eventlog, include_non_selfloops) {
	stop_eventlog(eventlog)
	case_classifier <- case_id(eventlog)
	size_of_selfloops <- size_of_selfloops_trace(eventlog) %>% select(-relative_trace_frequency)
	cases <- cases_light(eventlog) %>% select(-trace_id)

	r <- merge(cases, size_of_selfloops) %>% select(-trace) %>% tbl_df
	return(r)
}
