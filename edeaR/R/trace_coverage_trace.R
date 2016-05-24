
trace_coverage_trace <- function(eventlog) {
	stop_eventlog(eventlog)


	tr <- traces(eventlog)

	tr <- tr %>% select(trace, absolute_frequency, relative_frequency) %>% arrange(desc(absolute_frequency))
	tr$cum_sum <- cumsum(tr$relative_frequency)
	colnames(tr) <- c("trace","absolute","relative", "cum_sum")
	return(tr)
}
