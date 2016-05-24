
repetitions_activity <- function(eventlog) {
	stop_eventlog(eventlog)



	tr <- traces(eventlog)
	sl <- selfloops(eventlog)
	act <- activities(eventlog)
	ca <- cases_light(eventlog)

	colnames(sl)[colnames(sl)==activity_id(eventlog)] <- "event_classifier"
	colnames(act)[colnames(act)==activity_id(eventlog)] <- "event_classifier"

	ncases <- nrow(ca)

	colnames(eventlog)[colnames(eventlog)==activity_id(eventlog)] <- "event_classifier"
	colnames(eventlog)[colnames(eventlog)==case_id(eventlog)] <- "case_classifier"

	act_presence <- group_by(eventlog, event_classifier) %>%
		summarize(presence_in_cases = n_distinct(case_classifier))


	act <- merge(act, act_presence)
	sl <- merge(sl, ca)
	sl <- group_by(sl, event_classifier) %>%
		summarize(total_length = sum((length-1)))

	r <- merge(act, sl)
	r$number_of_repetitions <- r$absolute_frequency  - r$total_length - r$presence_in_cases
	r <- select(r, event_classifier, relative_frequency, number_of_repetitions) %>%
		mutate(relative = number_of_repetitions/ncases)

	colnames(r)[colnames(r)=="event_classifier"] <- activity_id(eventlog)
	colnames(r)[colnames(r)=="number_of_repetitions"] <- "absolute"
	r <- arrange(r, desc(relative_frequency))
	return(r)
}
