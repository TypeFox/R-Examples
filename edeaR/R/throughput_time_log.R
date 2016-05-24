


throughput_time_log <- function(eventlog,
								 units = "days") {
	stop_eventlog(eventlog)

	r <- throughput_time_case(eventlog, units = units)

	s <- summary(r$throughput_time)
	s <- c(s, St.Dev = sd(r$throughput_time))
	s <- c(s, IQR = s[5] - s[2])
	s <- c(s, tot = sum(r$throughput_time))
	names(s) <- c("min","q1","median","mean","q3","max","st_dev","iqr", "tot")

	s <- as.data.frame(s)
	s <- t(s)
	row.names(s) <- NULL
	return(s)
}
