

repetitions_trace <- function(eventlog) {

	stop_eventlog(eventlog)


	sl <- selfloops(eventlog)


	act2 <- select(activities(eventlog), -absolute_frequency)
	act <-act2[,1]
	colnames(act) <- "event_classifier"

	tra <- traces(eventlog)

	base <- merge(act, select(tra, trace))
	ca <- cases_light(eventlog)

	r <- merge(ca, eventlog)


	colnames(r)[colnames(r)==activity_id(eventlog)] <- "event_classifier"
	colnames(r)[colnames(r)==case_id(eventlog)] <- "case_classifier"
	colnames(r)[colnames(r)==activity_instance_id(eventlog)] <- "activity_instance_classifier"

	colnames(sl)[colnames(sl)==activity_id(eventlog)] <- "event_classifier"

	r <- group_by(r, trace, case_classifier,event_classifier) %>% summarize(activity_frequency = n_distinct(activity_instance_classifier))
	r <- group_by(r, trace, event_classifier) %>% summarize(activity_frequency = mean(activity_frequency))

	r <- merge(merge(base, r, all = T), act, all = T)
	r$activity_frequency[is.na(r$activity_frequency)] <- 0

	sl <- group_by(sl, trace, event_classifier) %>% summarize(total_length = sum((length-1)))

	r <- merge(r, sl, all.x = T)
	r$total_length[is.na(r$total_length)] <- 0
	r$number_of_repetitions <- pmax(0:0,(r$activity_frequency  - r$total_length - 1))

	r$instances <- r$activity_frequency - r$total_length

	r <- r %>% ungroup %>% group_by(trace) %>% summarize(absolute = sum(number_of_repetitions), denom = sum(instances))
	r$relative <- r$absolute / r$denom
	r <- r %>% select(-denom)
	r <- merge(select(tra, trace, relative_frequency),r)

	colnames(r)[colnames(r)=="relative_frequency"] <- "relative_trace_frequency"

	r <- r %>% arrange(desc(relative_trace_frequency)) %>% tbl_df()
	return(r)
}
