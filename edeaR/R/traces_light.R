traces_light <- function(eventlog){

	colnames(eventlog)[colnames(eventlog) == case_id(eventlog)] <- "case_classifier"
	colnames(eventlog)[colnames(eventlog) == activity_id(eventlog)] <- "event_classifier"
	colnames(eventlog)[colnames(eventlog) == timestamp(eventlog)] <- "timestamp_classifier"


	traces <- eventlog %>%
		group_by(case_classifier) %>%
		arrange(timestamp_classifier) %>%
		summarize(trace = paste(event_classifier, collapse = ",")) %>%
		group_by(trace) %>%
		summarize()

		return(traces)

}
