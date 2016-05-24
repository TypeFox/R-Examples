#' @title Filter: Time Period
#' @description Function to filter eventlog using a time period.
#' @param eventlog The event log to be used. An object of class
#' \code{eventlog}.
#' @param start_point Start timestamp of the time period. This should be a date object.
#' @param end_point End timestamp of the time period. This should be a data object.
#' @param filter_method Can be \code{contained, start, complete, intersecting} or \code{trim}.
#' \code{contained} keeps all the events related to cases contained in the time period.
#' \code{start} keeps all the events related to cases started in the time period.
#' \code{complete} keeps all the events related to cases complete in the time period.
#' \code{intersecting} keeps all the events related to cases in which at least one event started and/or ended in the time period.
#' \code{trim} keeps all the events which started and ended in the time frame.
#' @param reverse A logical parameter depicting whether the selection should be reversed.
#' @export

filter_time_period <- function(eventlog,
							   start_point,
							   end_point,
							   filter_method = "contained",
							   reverse = FALSE)
{
	stop_eventlog(eventlog)

	if(!(filter_method %in% c("contained","intersecting","start","complete","trim"))){
		stop(paste("Filter_method",
				   filter_method,
				   "non-valid. Should be one of: contained, intersecting, start, complete, trim.", sep = " "))
	}



	if(!("POSIXct" %in% class(start_point))) {
		stop("Start_point should be a date object.")
	}
	if(!("POSIXct" %in% class(end_point))) {
		stop("End_point should be a date object.")
	}

	c_sum <- cases(eventlog = eventlog)
	colnames(c_sum)[colnames(c_sum)==case_id(eventlog)] <- "case_classifier"
	colnames(eventlog)[colnames(eventlog)==case_id(eventlog)] <- "case_classifier"


	if(filter_method == "contained") {
		case_selection <- filter(c_sum, start_timestamp >= start_point, complete_timestamp <= end_point) %>%
			select(case_classifier)
	}
	else if(filter_method == "intersecting") {
		case_selection <- filter(c_sum, start_timestamp >= start_point & start_timestamp <= end_point |
								 	complete_timestamp >= start_point & complete_timestamp <= end_point |
								 	start_timestamp <= start_point & complete_timestamp >= end_point) %>%
			select(case_classifier)
	}
	else if(filter_method == "start") {
		case_selection <- filter(c_sum, start_timestamp >= start_point, start_timestamp <= end_point) %>%
			select(case_classifier)
	}
	else if(filter_method == "complete") {
		case_selection <- filter(c_sum, complete_timestamp >= start_point, complete_timestamp <= end_point) %>%
			select(case_classifier)
	}
	else if(filter_method == "trim") {



		colnames(eventlog)[colnames(eventlog) == case_id(eventlog)] <- "case_classifier"
		colnames(eventlog)[colnames(eventlog) == activity_id(eventlog)] <- "event_classifier"
		colnames(eventlog)[colnames(eventlog) == timestamp(eventlog)] <- "timestamp_classifier"
		colnames(eventlog)[colnames(eventlog) == activity_instance_id(eventlog)] <- "activity_instance_classifier"

		e <- eventlog %>%
			group_by(case_classifier, event_classifier, activity_instance_classifier) %>%
			summarize(start = min(timestamp_classifier), complete = max(timestamp_classifier))

		if(reverse == FALSE)
			f_eventlog <- filter(e, start >= start_point &
							   	complete <= end_point)
		else
			f_eventlog <- filter(e, !(start >= start_point &
							   	complete <= end_point))


		output <- filter(eventlog, activity_instance_classifier %in% f_eventlog$activity_instance_classifier)

		colnames(output)[colnames(output) == "case_classifier"] <- case_id(eventlog)
		colnames(output)[colnames(output) == "event_classifier"] <- activity_id(eventlog)
		colnames(output)[colnames(output) == "timestamp_classifier"] <- timestamp(eventlog)
		colnames(output)[colnames(output) == "activity_instance_classifier"] <- activity_instance_id(eventlog)

		return(eventlog(eventlog = output,
						activity_id = activity_id(eventlog),
						case_id = case_id(eventlog),
						timestamp =timestamp(eventlog),
						lifecycle_id = lifecycle_id(eventlog),
						activity_instance_id = activity_instance_id(eventlog)))
	}


	if(reverse == FALSE)
		f_eventlog <- filter(eventlog, (case_classifier %in% case_selection$case_classifier))
	else
		f_eventlog <- filter(eventlog, !(case_classifier %in% case_selection$case_classifier))

	colnames(f_eventlog)[colnames(f_eventlog) == "case_classifier"] <- case_id(eventlog)

	return(eventlog(eventlog = f_eventlog,
					activity_id = activity_id(eventlog),
					case_id = case_id(eventlog),
					timestamp =timestamp(eventlog),
					lifecycle_id = lifecycle_id(eventlog),
					activity_instance_id = activity_instance_id(eventlog)))
}
