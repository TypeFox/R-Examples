
end_activities_activity <- function(eventlog) {

	stop_eventlog(eventlog)

	tra <- cases_light(eventlog)

	for(i in 1:nrow(tra)){
		tra$last_event[i] <- strsplit(tra$trace[i], split = ",")[[1]][length(strsplit(tra$trace[i], split = ",")[[1]])]
	}

	ncases <- nrow(tra)

	r <- tra %>% group_by(last_event) %>% summarize(absolute = n()) %>% ungroup() %>% arrange(desc(absolute))
	r$relative <- r$absolute/ncases
	r$cum_sum <- cumsum(r$relative)
	colnames(r)[colnames(r) == "last_event"] <- activity_id(eventlog)
	return(r)

}
