#' @title Metric: Activity Presence
#'
#' @description Calculates for each activity type in what percentage of cases it is present.
#'
#' @param eventlog The event log to be used. An object of class
#' \code{eventlog}.
#'
#'
#' @seealso \code{\link{activity_type_frequency}}
#'
#' @examples
#'
#' data(example_log)
#' activity_presence(example_log)
#'
#'
#' @export activity_presence


activity_presence <- function(eventlog) {
	stop_eventlog(eventlog)

	event_classifier <- activity_id(eventlog)
	case_classifier <- case_id(eventlog)
	colnames(eventlog)[colnames(eventlog)==case_id(eventlog)] <- "case_classifier"

	r <- eventlog %>%
		group_by_(event_classifier, "case_classifier") %>%
		summarize() %>%
		summarize("absolute" = n_distinct(case_classifier)) %>%
		mutate(relative = absolute/n_cases(eventlog)) %>%
		arrange(desc(absolute))
	return(r)
}
