print_processing_times <- function(proc_times, start, end) {

	# Find total time
	times <- rep(NA, length(proc_times))
	for(le in 1:length(proc_times)) times[le] <- diff(proc_times[[le]])*1000
	total_time <- sum(times)

	#print(proc_times)
	
	cat('\nProcessing times\n')
	for(le in 1:length(proc_times)){
		cat(paste0('\t', names(proc_times)[le], ': ', round(times[le]), ' ms (', round((times[le]/total_time)*100), '%) \n'))
	}
	cat(paste0('\tTotal: ', round(total_time/1000, 2), ' s\n'))
	cat(paste0('\tTotal time: ', round(end-start, 2), ' s\n'))
	cat('\n')
}