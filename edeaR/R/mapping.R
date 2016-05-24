#' @title Mapping
#'
#' @description Prints the mapping of an event log object.
#'
#' @param eventlog The event log to be used. An object of class
#' \code{eventlog}.
#'
#' @export mapping

mapping <- function(eventlog) {
	stop_eventlog(eventlog)
	message(paste0("Case identifier:\t\t", case_id(eventlog)))
	message(paste0("Activity identifier:\t\t", activity_id(eventlog)))
	message(paste0("Activity instance identifier:\t", activity_instance_id(eventlog)))
	message(paste0("Timestamp:\t\t\t", timestamp(eventlog)))
	message(paste0("Lifecycle transition:\t\t", lifecycle_id(eventlog = eventlog)))
}
