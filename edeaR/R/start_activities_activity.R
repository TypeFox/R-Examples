

start_activities_activity <- function(eventlog) {

	stop_eventlog(eventlog)
	tra <- cases_light(eventlog)

	for(i in 1:nrow(tra)){
		tra$first_event[i] <- strsplit(tra$trace[i], split = ",")[[1]][1]
	}

	ncases <- nrow(tra)

	r <- tra %>% group_by(first_event) %>% summarize(absolute = n()) %>% ungroup() %>% arrange(desc(absolute))
	r$relative <- r$absolute/ncases
	r$cum_sum <- cumsum(r$relative)
	colnames(r)[colnames(r) == "first_event"] <- activity_id(eventlog)
	return(r)

}
