#' @title Metric: Number of traces
#'
#' @description Computes how many traces there are. The result is returned as absolute number as well as a relative number.
#' The relative number refers to the number of traces per 100 cases.
#'
#' @param eventlog The event log to be used. An object of class
#' \code{eventlog}.
#'
#'
#'
#' @export number_of_traces


number_of_traces <- function(eventlog) {

	stop_eventlog(eventlog)

	tr <- traces(eventlog)

	r<- data.frame( c(nrow(tr),sum(tr$absolute_frequency)/nrow(tr)))
	r <- as.data.frame(t(r))
	colnames(r) <- c("absolute","average_coverage")
	r <- tbl_df(r)
	return(r)
}
