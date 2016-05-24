

repetitions_log <- function(eventlog) {

	stop_eventlog(eventlog)

	rep <- repetitions_case(eventlog)

	s <- summary(rep$absolute)
	s <- c(s, St.Dev = sd(rep$absolute))
	s <- c(s, IQR = s[5] - s[2])
	s <- c(s, tot = sum(rep$absolute))
	names(s) <- c("min","q1","median","mean","q3","max","st_dev","iqr", "tot")

	s <- as.data.frame(s)
	s <- t(s)
	row.names(s) <- NULL
	return(s)

}
