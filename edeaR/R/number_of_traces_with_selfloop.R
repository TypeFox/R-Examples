#' @title Metric: Number of Traces with Selfloop
#'
#' @description Returns the number of traces in which one or more selfloops occur, both in absolute and relative numbers.
#'
#'
#' @param eventlog The event log to be used. An object of class
#' \code{eventlog}.
#'
#'
#' @export number_of_traces_with_selfloop

number_of_traces_with_selfloop <- function(eventlog) {

	stop_eventlog(eventlog)
	r <- number_of_selfloops(eventlog,	   level_of_analysis = "trace")
	ntraces <- nrow(r)
	r <- r %>% filter(absolute > 0)

	a <- nrow(r)
	r <- data.frame(c(a, a/ntraces))
	r <- as.data.frame(t(r))
	colnames(r) <- c("absolute","relative")
	row.names(r) <- NULL
	r <- tbl_df(r)
	return(r)
}
