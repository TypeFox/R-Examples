
trace_coverage_log <- function(eventlog,
						   threshold = NULL) {
	stop_eventlog(eventlog)
	if(is.null(threshold)) {
		warning("Threshold defaulted to 0.8. Use `threshold = x` to adjust this")
		threshold = 0.8
	}

	tra <- traces(eventlog)

	tr <- tra %>% arrange(desc(relative_frequency)) %>% group_by(relative_frequency) %>% summarize(s = sum(relative_frequency)) %>% arrange(desc(relative_frequency))
	tr$c <- cumsum(tr$s)
	if(tr$c[1] >= threshold){
		tr <- tr[1,]
	}
	else if(threshold == 1)
		tr <- tr[nrow(tr),]
	else {
		stop = FALSE
		for(i in 2:nrow(tr)){
			if(!stop && tr$c[i-1] <= threshold && tr$c[i] >= threshold){
				tr <- tr[(i-1):i,]
				stop = TRUE
			}
		}
	}
	tr$cnt <- 0:0
	for(i in 1:nrow(tr))
		tr$cnt[i] <- nrow(filter(tra, relative_frequency >= tr$relative_frequency[i] ))
	tr <- tr %>% select(cnt,c)
	colnames(tr) <- c("number_of_traces","coverage")
	tr <- tbl_df(tr)
	return(tr)

}
