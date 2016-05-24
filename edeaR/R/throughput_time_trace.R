

throughput_time_trace <- function(eventlog, units = "days") {
	stop_eventlog(eventlog)

	case_classifier <- case_id(eventlog)

	t <- throughput_time_case(eventlog, units = units)
	c <- cases_light(eventlog)

	r <- merge(t, c, by = case_classifier) %>%
		group_by(trace) %>%
		summarize(relative_trace_frequency = n(),
				  min = min(throughput_time),
				  q1 = quantile(throughput_time, probs = c(0.25)),
				  median = median(throughput_time),
				  mean = mean(throughput_time),
				  q3 = quantile(throughput_time, probs = c(0.75)),
				  max = max(throughput_time),
				  st_dev = sd(throughput_time),
				  iqr = quantile(throughput_time, probs = c(0.75)) - quantile(throughput_time, probs = c(0.25)),
				  tot = sum(throughput_time)) %>%
		mutate(relative_trace_frequency = relative_trace_frequency/sum(relative_trace_frequency)) %>%
		arrange(desc(relative_trace_frequency))

	return(r)


}
