

processing_time_trace <- function(eventlog,
								  units = "days") {
	stop_eventlog(eventlog)

	colnames(eventlog)[colnames(eventlog) == timestamp(eventlog)] <- "timestamp_classifier"
	case_classifier <- case_id(eventlog)
	event_classifier <- activity_id(eventlog)
	activity_instance_classifier <- activity_instance_id(eventlog)

	e <- eventlog %>%
		group_by_(case_classifier, activity_instance_classifier) %>%
		summarize(s = min(timestamp_classifier), e = max(timestamp_classifier)) %>%
		mutate(dur = as.double(e - s, units = units)) %>%
		summarize(tot_dur = sum(dur))

	ca <- cases_light(eventlog)

	r <- merge(ca,e, by = case_classifier) %>%
		select(trace, tot_dur) %>%
		group_by(trace) %>%
		summarize(relative_frequency = n(),
				  min = min(tot_dur),
				  q1 = quantile(tot_dur, probs = c(0.25)),
				  median = median(tot_dur),
				  mean = mean(tot_dur),
				  q3 = quantile(tot_dur, probs = c(0.75)),
				  max = max(tot_dur),
				  st_dev = sd(tot_dur),
				  iqr = quantile(tot_dur, probs = c(0.75)) - quantile(tot_dur,probs = c(0.25)),
				  tot = sum(tot_dur)) %>%
		mutate(relative_frequency = relative_frequency/sum(relative_frequency)) %>%
		arrange(desc(relative_frequency))

	return(r)

}
