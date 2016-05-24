


filter_endpoints_percentile <- function(eventlog,
							 percentile_cut_off,
							 reverse = F) {


	stop_eventlog(eventlog)

	s_act <- start_activities_activity(eventlog)
	e_act <- end_activities_activity(eventlog)

	start_activities <- filter(s_act, cum_sum <= percentile_cut_off)
	end_activities <- filter(e_act, cum_sum <= percentile_cut_off)

	colnames(start_activities)[colnames(start_activities) == activity_id(eventlog)] <- "event_Classifier"
	colnames(end_activities)[colnames(end_activities) == activity_id(eventlog)] <- "event_Classifier"

	c_sum <- cases(eventlog = eventlog)
	colnames(c_sum)[colnames(c_sum)==case_id(eventlog)] <- "case_classifier"
	colnames(eventlog)[colnames(eventlog)==case_id(eventlog)] <- "case_classifier"


	case_selection <- filter(c_sum,
							 first_event %in% start_activities$event_classifier,
							 last_event %in% end_activities$event_classifier)$case_classifier

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
