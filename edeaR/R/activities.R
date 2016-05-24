
#' @title Activities
#'
#' @description Returns a \code{tbl_df}  containing a list of all activity types in the event log, with there absolute and relative frequency
#'
#' @param eventlog The event log to be used. An object of class
#' \code{eventlog}.
#'
#'
#' @seealso \code{\link{activity_id}},\code{\link{activity_instance_id}}, \code{\link{eventlog}}
#'
#' @examples
#'
#' data(example_log)
#' activities(example_log)
#'
#' @export activities


activities <- function(eventlog) {
	stop_eventlog(eventlog)
	colnames(eventlog)[colnames(eventlog) == activity_id(eventlog)] <- "event_classifier"
	colnames(eventlog)[colnames(eventlog) == activity_instance_id(eventlog)] <- "activity_instance_classifier"

	output <- eventlog %>%
		group_by_("event_classifier") %>%
		summarize("absolute_frequency" = n_distinct(activity_instance_classifier)) %>%
		arrange(absolute_frequency)
	output$relative_frequency <- output$absolute_frequency/sum(output$absolute_frequency)

	colnames(output)[colnames(output)== "event_classifier"] <- activity_id(eventlog)
	return(output)
}
