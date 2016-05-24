
filter_trace_length_percentile <- function(eventlog,
								percentile_cut_off,
								reverse = F)
{

	colnames(eventlog)[colnames(eventlog)==case_id(eventlog)] <- "case_classifier"
	colnames(eventlog)[colnames(eventlog)==activity_instance_id(eventlog)] <- "activity_instance_classifier"

	case_lengths <- group_by(eventlog, case_classifier) %>% summarize(length = n_distinct(activity_instance_classifier))
	case_lengths_grouped <- group_by(case_lengths, length) %>% summarize(freq = n()) %>% arrange(length)
	case_lengths_grouped$perc <- cumsum(case_lengths_grouped$freq)/sum(case_lengths_grouped$freq)
	case_lengths <- merge(case_lengths, case_lengths_grouped)

	case_selection <- filter(case_lengths, perc <= percentile_cut_off)$case_classifier

		if(reverse == FALSE)
		f_eventlog <- filter(eventlog, case_classifier %in% case_selection)
	else
		f_eventlog <- filter(eventlog, !(case_classifier %in% case_selection))

	colnames(f_eventlog)[colnames(f_eventlog)=="case_classifier"] <- case_id(eventlog)
	colnames(f_eventlog)[colnames(f_eventlog)=="activity_instance_classifier"] <- activity_instance_id(eventlog)

	output <- eventlog(f_eventlog,
					   activity_id = activity_id(eventlog),
					   case_id = case_id(eventlog),
					   timestamp =timestamp(eventlog),
					   lifecycle_id = lifecycle_id(eventlog),
					   activity_instance_id = activity_instance_id(eventlog))

	return(output)

}
