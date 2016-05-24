
repetitions_case <- function(eventlog) {
	stop_eventlog(eventlog)

	case_classifier <- case_id(eventlog)
	repetitions <- repetitions_trace(eventlog) %>% select(trace, absolute, relative)
	cases <- cases_light(eventlog) %>% select(-trace_id)

	r <- merge(repetitions, cases) %>%
		select_(case_classifier, "absolute","relative") %>%
		arrange(-absolute) %>%
		tbl_df

	return(r)
}
