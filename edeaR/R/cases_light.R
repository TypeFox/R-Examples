
cases_light <- function(eventlog){
	if(!("eventlog" %in% class(eventlog)))
		stop("Function only applicable for eventlog object")

	traces(eventlog, output_traces = FALSE, output_cases = TRUE)
}
