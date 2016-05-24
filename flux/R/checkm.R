checkm <- function(modelled, measured, t.unit = NULL){
	if(is.null(t.unit)){
		mots <- modelled$timestamp
		mets <- measured$timestamp
		INDEX <- apply(outer(mots, mets, "-"), 2, function(x) which.min(abs(x)))
		SECIND <- abs(difftime(mots[INDEX], mets, "secs")) < 3600
		INDEX <- INDEX[SECIND]
		tcheck <- data.frame(modelled[INDEX,], measured[SECIND,])
	}
	else{
		# make timestamp a character vector (for merging)
		modelled$ts <- as.character(modelled$timestamp)
		# round timestamp to next half hour
		measured$timestamp.round <- round(as.POSIXlt(measured$timestamp), t.unit)
		# make timestamp a character vector for merging
		measured$ts <- as.character(measured$timestamp.round)
		# merge and keep only data rows that exist in both tables
		tcheck <- merge(measured, modelled, by="ts")
	}
	# return merged data
	return(tcheck)
}