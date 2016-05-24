activity_type_frequency_case <- function(eventlog) {

	stop_eventlog(eventlog)

	case_classifier <- case_id(eventlog)
	event_classifier <- activity_id(eventlog)

	colnames(eventlog)[colnames(eventlog) == activity_instance_id(eventlog)] <- "activity_instance_classifier"
	ca <- cases_light(eventlog)
	e <- eventlog

	s <- e %>%
		group_by_(case_classifier, event_classifier) %>%
		summarize(Freq = n_distinct(activity_instance_classifier)) %>%
		summarize(min = min(Freq),
				  q1 =quantile(Freq, probs=0.25),
				  mean = mean(Freq),
				  median = median(Freq) * 1.0,
				  q3 =quantile(Freq, probs=0.75),
				  max = max(Freq),
				  st_dev = sd(Freq),
				  iqr = quantile(Freq, probs=0.75) - quantile(Freq, probs=0.25)
		)



	return(s)


}
