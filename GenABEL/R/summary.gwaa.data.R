"summary.gwaa.data" <- 
		function(object,...) {
	#ret <- list(phdata = summary(object@phdata),
	#	   gtdata = summary(object@gtdata))
	#ret
	return(summary(object@gtdata))
}
