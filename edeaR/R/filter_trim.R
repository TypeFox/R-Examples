#' @title Filter: Trim cases
#'
#' @description Trim all cases from the first event of a set of start activities to the last event of a set of end activities.
#' Traces that don't have at least one event of both sets are discarded.
#'
#' @param eventlog The event log to be used. An object of class
#' \code{eventlog}.
#'
#' @param start_activities Start activities used for trimming. If not provided, the start of the cases is not trimmed.
#'
#' @param end_activities End activities used for trimming. If not provided, the end of the cases or not trimmed.
#'
#' @param reverse A logical parameter depicting whether the selection should be reversed.
#'
#' @export filter_trim


filter_trim <- function(eventlog,
						start_activities = NULL,
						end_activities = NULL,
						reverse = F) {

	stop_eventlog(eventlog)

	if(is.null(start_activities) & is.null(end_activities))
		stop("At least on start or end activity should be provided")



	if(is.null(start_activities))
		start_activities <- unique(e$event_classifier)
	if(is.null(end_activities))
		end_activities <- unique(e$event_classifier)



	colnames(eventlog)[colnames(eventlog) == case_id(eventlog)] <- "case_classifier"
	colnames(eventlog)[colnames(eventlog) == activity_id(eventlog)] <- "event_classifier"
	colnames(eventlog)[colnames(eventlog) == timestamp(eventlog)] <- "timestamp_classifier"
	colnames(eventlog)[colnames(eventlog) == activity_instance_id(eventlog)] <- "activity_instance_classifier"


	ranked_eventlog <- group_by(eventlog, case_classifier, activity_instance_classifier, event_classifier) %>%
		summarize(min_timestamp = min(timestamp_classifier)) %>%
		group_by(case_classifier) %>%
		mutate(r = dense_rank(min_timestamp))
	min_ranks <- ranked_eventlog %>% filter(event_classifier %in% start_activities) %>% summarize(min_rank = min(r))
	max_ranks <- ranked_eventlog %>% filter(event_classifier %in% end_activities) %>% summarize(max_rank = max(r))
	ranked_eventlog <- merge(merge(ranked_eventlog, min_ranks), max_ranks)
	ranked_eventlog$included <- ranked_eventlog$r >= ranked_eventlog$min_rank & ranked_eventlog$r <= ranked_eventlog$max_rank


	if(reverse == F)
		ranked_eventlog <- ranked_eventlog[ranked_eventlog$r >= ranked_eventlog$min_rank & ranked_eventlog$r <= ranked_eventlog$max_rank,]
	else
		ranked_eventlog <- ranked_eventlog[!(ranked_eventlog$r >= ranked_eventlog$min_rank & ranked_eventlog$r <= ranked_eventlog$max_rank),]
	ranked_eventlog <- select(ranked_eventlog, -r, -min_rank, -max_rank)



	f_eventlog <- eventlog[eventlog$activity_instance_classifier %in% ranked_eventlog$activity_instance_classifier,]

	colnames(f_eventlog)[colnames(f_eventlog) == "case_classifier"] <- case_id(eventlog)
	colnames(f_eventlog)[colnames(f_eventlog) == "event_classifier"] <- activity_id(eventlog)
	colnames(f_eventlog)[colnames(f_eventlog) == "timestamp_classifier"] <- timestamp(eventlog)
	colnames(f_eventlog)[colnames(f_eventlog) == "activity_instance_classifier"] <- activity_instance_id(eventlog)


	eventlog <- eventlog(eventlog = f_eventlog,
						 activity_id = activity_id(eventlog),
						 case_id = case_id(eventlog),
						 activity_instance_id = activity_instance_id(eventlog),
						 lifecycle_id = lifecycle_id(eventlog),
						 timestamp = timestamp(eventlog))

	return(eventlog)
}
