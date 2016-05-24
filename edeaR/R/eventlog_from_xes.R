#' @title eventlog_from_xes
#' @description Extracts eventlog from a xes-file.
#' @param xesfile Reference to a .xes file, conforming to the xes-standard.
#' @seealso \url{http://www.xes-standard.org/}
#' @export eventlog_from_xes

eventlog_from_xes <- function(xesfile = file.choose()){
	log <- csv_from_xes(xesfile)
	colnames(log) <- gsub(":",".", colnames(log))
	log <- arrange(log, case_concept.name, event_concept.name, event_time.timestamp)

	log$activity_instance[1] <- 1
	for(i in 2:nrow(log))

		if(log$event_lifecycle.transition[i-1] %in% c("autoskip","manualskip","complete","withdraw","abort_activity","abort_case"))
			log$activity_instance[i] <- log$activity_instance[i-1] + 1
		else
			log$activity_instance[i] <- log$activity_instance[i-1]

	elog <- eventlog(eventlog = log,
					 case_id = "case_concept.name",
					 activity_id = "event_concept.name",
					 activity_instance_id = "activity_instance",
					 lifecycle_id = "event_lifecycle.transition",
					 timestamp = "event_time.timestamp")
	return(elog)
}
