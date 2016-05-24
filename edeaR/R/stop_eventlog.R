stop_eventlog <- function(eventlog)
	if(!("eventlog" %in% class(eventlog)))
		stop("Function only applicable for class eventlog")
