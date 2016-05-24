

trace_length_trace <- function(eventlog,
							   threshold = 0.8) {

	stop_eventlog(eventlog)

	tra <- traces(eventlog)

	tra_length <- group_by(tra, trace_id,
								   trace,
								   absolute_frequency) %>%
		summarize(number_of_events = length(strsplit(trace, split = ",")[[1]]))

	tra_freq <- group_by(tra, absolute_frequency) %>%
		summarize(s = sum(relative_frequency)) %>%
		arrange(desc(absolute_frequency))

	tra_freq$c <- cumsum(tra_freq$s)

	if(tra_freq$c[1] >= threshold){
		tr <- tra_freq[1,]
	}
	else if(threshold == 1) {
		tr <- tra_freq[nrow(tra_freq),]
	}
	else {
		stop = FALSE
		for(i in 2:nrow(tra_freq)){
			if(!stop && tra_freq$c[i-1] <= threshold && tra_freq$c[i] >= threshold){
				tr <- tra_freq[(i-1):i,]
				stop = TRUE
			}
		}
	}


	for(i in 1:nrow(tr)){
		f <- tra_length %>% filter(absolute_frequency >= tr$absolute_frequency[i])
		avg_length <- sum(f$absolute_frequency * f$number_of_events) / sum(f$absolute_frequency)
		tra_length[,paste("relative_to_top_", round(tr$c[i]*100, 2), sep = "")] <- tra_length$number_of_events/avg_length
	}

	tra_length <- merge(select(tra, trace, relative_frequency), tra_length)
	tra_length <- tra_length %>%  arrange(desc(absolute_frequency)) %>% select(-trace_id, -absolute_frequency)

	colnames(tra_length)[colnames(tra_length) == "number_of_events"] <- "absolute"
	colnames(tra_length)[colnames(tra_length) == "relative_frequency"] <- "relative_trace_frequency"

	tra_length <- tbl_df(tra_length)
	return(tra_length)

}
