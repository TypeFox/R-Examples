
processing_time_activity <- function(eventlog,
									 units = "days") {
	stop_eventlog(eventlog)

	colnames(eventlog)[colnames(eventlog) == timestamp(eventlog)] <- "timestamp_classifier"
	event_classifier <- activity_id(eventlog)
	activity_instance_classifier <- activity_instance_id(eventlog)

	r <- eventlog %>%
		group_by_(event_classifier, activity_instance_classifier) %>%
		summarize(s = min(timestamp_classifier), e = max(timestamp_classifier)) %>%
		mutate(dur = as.double(e - s, units = units)) %>%
		summarize(relative_frequency = n(),
				  min = min(dur),
				  q1 = quantile(dur, probs = c(0.25)),
				  median = median(dur),
				  mean = mean(dur),
				  q3 = quantile(dur, probs = c(0.75)),
				  max = max(dur),
				  st_dev = sd(dur),
				  iqr = quantile(dur, probs = c(0.75)) - quantile(dur,probs = c(0.25)),
				  tot = sum(dur)) %>%
		mutate(relative_frequency = relative_frequency/sum(relative_frequency)) %>%
		arrange(desc(relative_frequency))

	return(r)
}
