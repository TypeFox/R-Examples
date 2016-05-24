activity_type_frequency_trace <- function(eventlog) {

	stop_eventlog(eventlog)


	case_classifier <- case_id(eventlog)
	event_classifier <- activity_id(eventlog)
	colnames(eventlog)[colnames(eventlog) == activity_instance_id(eventlog)] <- "activity_instance_classifier"
	ca <- cases_light(eventlog)
	e <- eventlog

	s <- e %>%
		group_by_(case_classifier, event_classifier) %>%
		summarize(Freq = n_distinct(activity_instance_classifier)) %>%
		group_by_(case_classifier) %>%
		summarize(min = min(Freq),
				  q1 =quantile(Freq, probs=0.25),
				  mean = mean(Freq),
				  median = median(Freq) * 1.0,
				  q3 =quantile(Freq, probs=0.75),
				  max = max(Freq),
				  st_dev = sd(Freq),
				  iqr = quantile(Freq, probs=0.75) - quantile(Freq, probs=0.25)
		)

	r <- merge(s, ca)

	r <- r %>%
		group_by(trace, min, q1, mean, median, q3, max, st_dev) %>%
		summarize(relative_trace_frequency = n()) %>%
		ungroup %>%
		select(trace, relative_trace_frequency, min, q1, mean, median, q3, max, st_dev) %>%
		mutate(iqr = q3 - q1) %>%
		arrange(desc(relative_trace_frequency))


	r$relative_trace_frequency <- r$relative_trace_frequency/(sum(r$relative_trace_frequency))
	return(r)


}
