


filter_trace_frequency_threshold <- function(eventlog,
											  lower_threshold = NULL,
											  upper_threshold = NULL,
											  reverse = F){

	if(is.null(lower_threshold) & is.null(upper_threshold)){
		stop("Upper threshold or lower threshold must be defined")
	}
	if(is.null(lower_threshold))
		lower_threshold <- -Inf
	if(is.null(upper_threshold))
		upper_threshold <- Inf


	if(reverse == F)
		case_selection <- merge(cases_light(eventlog),
								trace_frequency(eventlog) %>%
									filter(absolute >= lower_threshold,
										   absolute <= upper_threshold))
	else
		case_selection <- merge(cases_light(eventlog),
								trace_frequency(eventlog) %>%
									filter(absolute < lower_threshold |
										   absolute > upper_threshold))


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
