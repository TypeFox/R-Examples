size_of_selfloops_log <- function(eventlog, include_non_selfloops = FALSE) {


	stop_eventlog(eventlog)

	sl <- selfloops(eventlog)
	ca <- cases_light(eventlog)
	sl <- merge(sl, ca)


	selfloops_lengths <- sl$length
	if(!include_non_selfloops){
		selfloops_lengths <- selfloops_lengths - 1

		if(length(selfloops_lengths) > 0){
			r <- summary(selfloops_lengths)
			r <- c(r, St.Dev = sd(selfloops_lengths))
			r <- c(r, IQR = r[5] - r[2])
		}
		else {
			r <- rep(NA, 8)
		}
	}
	else {
		number_of_events <- nrow(eventlog)
		number_of_zero_selfloops <- number_of_events - sum(selfloops_lengths)

		selfloops_lengths <- c(selfloops_lengths, rep(1, number_of_zero_selfloops))
		selfloops_lengths <- selfloops_lengths - 1

		r <- summary(selfloops_lengths)
		r <- c(r, St.Dev = sd(selfloops_lengths))
		r <- c(r, IQR = r[5] - r[2])
	}


	names(r) <- c("min","q1","median","mean","q3","max","st_dev","iqr")
	r <- data.frame(r)
	r <- as.data.frame(t(r))

	row.names(r) <- NULL

	r <- tbl_df(r)
	return(r)
}
