
filter_trace_length_threshold <- function(eventlog,
									  lower_threshold = NULL,
									  upper_threshold = NULL,
									  reverse = F)
{
	if(is.null(lower_threshold) & is.null(upper_threshold)){
		stop("Upper threshold or lower threshold must be defined")
	}

	if(is.null(lower_threshold))
		lower_threshold <- -Inf
	if(is.null(upper_threshold))
		upper_threshold <- Inf

	colnames(eventlog)[colnames(eventlog)==case_id(eventlog)] <- "case_classifier"
	colnames(eventlog)[colnames(eventlog)==activity_instance_id(eventlog)] <- "activity_instance_classifier"

	case_lengths <- group_by(eventlog, case_classifier) %>% summarize(length = n_distinct(activity_instance_classifier))

	case_selection <- filter(case_lengths, length >= lower_threshold, length <= upper_threshold)$case_classifier
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
