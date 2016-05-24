

end_activities_log <- function(eventlog) {

	stop_eventlog(eventlog)

	r <- end_activities_activity(eventlog)

	n_activities <- nrow(activities(eventlog))

	r<- data.frame(c(nrow(r),nrow(r)/n_activities))
	r <- as.data.frame(t(r))
	colnames(r) <- c("absolute","relative")
	row.names(r) <- NULL
	r <- tbl_df(r)
	return(r)
}
