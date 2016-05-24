



generate_pattern_dummies <- function(patterns,
									 eventlog,
									 interleavings_allowed = TRUE) {
	stop_eventlog(eventlog)

	cases <- cases_light(eventlog)
	colnames(cases)[colnames(cases)==case_id(eventlog)] <- "case_classifier"
	variants <- traces(eventlog)
	variants <- arrange(variants, trace_id)

	dummies <- list()
	empty_vector <- c(seq(from = 0, to = 0, length.out = nrow(patterns)))

	for(i in 1:nrow(variants)){
		dummies[[i]] <-  empty_vector
	}
	for(i in 1:nrow(variants)) {
		trace <- strsplit(toString(variants$trace[i]),split = ",")[[1]]
		for(j in 1:nrow(patterns)) {
			pattern <- strsplit(toString(patterns$pattern[j]),split =",")[[1]]
			if(contains_sequence(pattern, trace, interleavings_allowed = interleavings_allowed))
				dummies[[i]][j] <- 1
		}
	}

	patterndummies_per_case <- data.frame(matrix(nrow = nrow(cases), ncol = nrow(patterns),0:0))
	colnames(patterndummies_per_case) <- paste("X", c(1:nrow(patterns)), sep = "")
	patterndummies_per_case$case_classifier <- cases$case_classifier
	for(i in 1:nrow(patterndummies_per_case)){

		pointer <- cases$trace_id[i]
		for(j in 1:nrow(patterns)){
			patterndummies_per_case[i,j] <- dummies[[pointer]][j]
		}
	}
	output<- merge(patterndummies_per_case,cases, by = colnames(cases)[1])
	colnames(output)[colnames(output)=="case_classifier"] <- case_id(eventlog)
	return(output)
}
