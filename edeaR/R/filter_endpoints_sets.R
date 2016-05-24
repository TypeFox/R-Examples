

filter_endpoints_sets <- function(eventlog,
							 start_activities = NULL,
							 end_activities = NULL,
							 reverse = F) {

	stop_eventlog(eventlog)

	if(is.null(start_activities) & is.null(end_activities)){
		stop("At least on start or end activity for filtering should be provided")
	}

	if(is.null(start_activities))
		start_activities <- unique(eventlog$event_classifier)
	if(is.null(end_activities))
		end_activities <- unique(eventlog$event_classifier)

	c_sum <- cases(eventlog = eventlog)
	colnames(c_sum)[colnames(c_sum)==case_id(eventlog)] <- "case_classifier"
	colnames(eventlog)[colnames(eventlog)==case_id(eventlog)] <- "case_classifier"

	case_selection <- filter(c_sum, first_event %in% start_activities, last_event %in% end_activities)$case_classifier

	if(reverse == FALSE)
		f_eventlog <- filter(eventlog, case_classifier %in% case_selection)
	else
		f_eventlog <- filter(eventlog, !(case_classifier %in% case_selection))

	colnames(f_eventlog)[colnames(f_eventlog)=="case_classifier"] <- case_id(eventlog)
	colnames(f_eventlog)[colnames(f_eventlog)=="event_classifier"] <- activity_id(eventlog)


	output <- eventlog(eventlog = f_eventlog,
					   activity_id = activity_id(eventlog),
					   case_id = case_id(eventlog),
					   timestamp =timestamp(eventlog),
					   lifecycle_id = lifecycle_id(eventlog),
					   activity_instance_id = activity_instance_id(eventlog))

	return(output)

}
