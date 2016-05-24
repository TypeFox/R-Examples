
filter_trace_frequency_percentile <- function(eventlog,
								  percentile_cut_off,
								  reverse = F){




	if(reverse == F)
		case_selection <- merge(cases_light(eventlog),
								trace_frequency(eventlog) %>%
									filter(cum_sum <= percentile_cut_off))
	else
		case_selection <- merge(cases_light(eventlog),
								trace_frequency(eventlog) %>%
									filter(cum_sum > percentile_cut_off))



	colnames(eventlog)[colnames(eventlog)==case_id(eventlog)] <- "case_classifier"
	colnames(case_selection)[colnames(case_selection) == case_id(eventlog)] <- "case_classifier"


	output <- filter(eventlog,  case_classifier %in% case_selection$case_classifier)
	colnames(output)[colnames(output) == "case_classifier"] <- case_id(eventlog)

	output <- eventlog(output,
					   activity_id = activity_id(eventlog),
					   case_id = case_id(eventlog),
					   timestamp =timestamp(eventlog),
					   lifecycle_id = lifecycle_id(eventlog),
					   activity_instance_id = activity_instance_id(eventlog))

	return(output)

}
