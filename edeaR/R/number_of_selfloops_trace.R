

number_of_selfloops_trace <- function(eventlog) {

	stop_eventlog(eventlog)

	tr <- traces(eventlog)
	sl <- selfloops(eventlog)

	output <- merge(tr, sl)
	output <- output %>% group_by(trace, relative_frequency) %>% summarize(absolute = n(), ev_in = sum(length),n_ev = length(strsplit(trace, split = ",")[[1]])) %>% ungroup()

	output$ni <- output$n_ev - output$ev_in + output$absolute

	output$relative <- output$absolute/output$ni

	output <- merge(output, select(tr, trace, relative_frequency), all.y = T)

	output[is.na(output$absolute), "absolute"] <- 0
	output[is.na(output$relative), "relative"] <- 0

	output <- output %>% select(trace, relative_frequency, absolute, relative) %>% arrange( desc(relative_frequency))

	colnames(output)[colnames(output)=="relative_frequency"] <- "relative_trace_frequency"
	output <- tbl_df(output)
	return(output)
}



