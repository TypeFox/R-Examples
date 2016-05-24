

selfloops <- function(eventlog) {

	stop_eventlog(eventlog)

	tr <- traces(eventlog)
	output <- data.frame(trace = character(0), length = numeric(0), Activity = character(0) , stringsAsFactors = FALSE)
	r <- 1
	for(i in 1:nrow(tr)) {

		act_seq <- strsplit(tr$trace[i], split = ",")[[1]]
		if(length(act_seq) > 1){
			current_act <- act_seq[1]
			length <- 1
			for(j in 2:length(act_seq)){
				if(current_act == act_seq[j])
					length <- length + 1
				else {
					if(length > 1){
						output <- bind_rows(output, data.frame(t(c(trace = tr$trace[i],length = 9999, Activity = current_act))))
						output$length[r] <- length
						r <- r + 1
					}
					length <- 1
					current_act <- act_seq[j]
				}
			}
			if(length > 1){
				output <- bind_rows(output, data.frame(t(c(trace = tr$trace[i],length = 9999, Activity = current_act))))
				output$length[r] <- length
				r <- r + 1
			}
		}
	}
	output <- tbl_df(output)
	colnames(output)[colnames(output) == "Activity"] <- activity_id(eventlog)
	return(output)
}
