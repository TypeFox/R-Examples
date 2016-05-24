

number_of_selfloops_activity <- function(eventlog) {

	stop_eventlog(eventlog)

	sl <- selfloops(eventlog)
	tr <- traces(eventlog)
	act <- activities(eventlog) %>% rename(act_freq = absolute_frequency) %>% select(1:2)
	r <- merge(sl, tr, "trace")

	colnames(r)[colnames(r) == activity_id(eventlog)] <- "event_classifier"

	r <- group_by(r, event_classifier) %>% summarize(absolute = sum(absolute_frequency), total_length = sum(absolute_frequency*(length-1)))

	colnames(r)[colnames(r) == "event_classifier"] <- activity_id(eventlog)

	r <- merge(r, act)
	r <- r %>% mutate(relative = absolute/(act_freq - total_length)) %>% select(-total_length, -act_freq)
	r <- arrange(r, desc(absolute))
	return(r)

}
