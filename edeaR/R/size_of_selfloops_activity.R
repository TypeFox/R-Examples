


size_of_selfloops_activity <- function(eventlog, include_non_selfloops = FALSE) {

	stop_eventlog(eventlog)


	act <- activities(eventlog)
	sl <- selfloops(eventlog)
	colnames(act)[colnames(act) == activity_id(eventlog)] <- "event_classifier"
	colnames(sl)[colnames(sl) ==   activity_id(eventlog)] <- "event_classifier"

	for(i in 1:nrow(act)){
		selfloops_lengths <- filter(sl, event_classifier == act$event_classifier[i])$length

		if(!include_non_selfloops){
			selfloops_lengths <- selfloops_lengths - 1

			if(length(selfloops_lengths) > 0){
				act$min[i] <- min(selfloops_lengths)
				act$q1[i] <- quantile(selfloops_lengths, probs = 0.25)
				act$mean[i] <- mean(selfloops_lengths)
				act$median[i] <- median(selfloops_lengths)
				act$q3[i] <- quantile(selfloops_lengths, probs = 0.75)
				act$max[i] <- max(selfloops_lengths)
				act$st_dev[i] <- sd(selfloops_lengths)
				act$iqr[i] <- act$q3[i] - act$q1[i]
			}
			else {
				act$min[i] <- NA
				act$q1[i] <- NA
				act$mean[i] <- NA
				act$median[i] <- NA
				act$q3[i] <- NA
				act$max[i] <- NA
				act$st_dev[i] <- NA
				act$iqr[i] <- NA}
		}
		else {

			number_of_zero_selfloops <- act$absolute_frequency[i] - sum(selfloops_lengths)
			selfloops_lengths <- c(selfloops_lengths, rep(1, number_of_zero_selfloops))
			selfloops_lengths <- selfloops_lengths - 1
			act$min[i] <- min(selfloops_lengths)
			act$q1[i] <- quantile(selfloops_lengths, probs = 0.25)
			act$mean[i] <- mean(selfloops_lengths)
			act$median[i] <- median(selfloops_lengths)
			act$q3[i] <- quantile(selfloops_lengths, probs = 0.75)
			act$max[i] <- max(selfloops_lengths)
			act$st_dev[i] <- sd(selfloops_lengths)
			act$iqr[i] <- act$q3[i] - act$q1[i]
		}
	}

	act <- act %>% select(-absolute_frequency) %>% arrange(desc(relative_frequency))
	colnames(act)[colnames(act)== "relative_frequency"] <- "relative_activity_frequency"
	colnames(act)[colnames(act)=="event_classifier"] <- activity_id(eventlog)
	return(act)
}
