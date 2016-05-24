pcitMemoryRequirement <- function(m, units=c("MB", "bytes", "KB", "GB", "TB"), nCopies=3) {
	double.bytes <- 8
	units <- match.arg(units)
	
	if (class(m)=="numeric" | class(m)=="integer") {
		nNodes <- m
	} else {
		nNodes <- try(nrow(m), silent=TRUE)
		if (class(nNodes) == "try-error" | is.null(nNodes)) {
			cat("ERROR: argument 'm' must be a numeric OR an object on which nrow() can be performed.\n\n", geterrmessage())
			return(FALSE)
		}
	}
	
	switch(units,
			"bytes" = { denominator <- 1024^0 },
			"KB" = { denominator <- 1024^1 },
			"MB" = { denominator <- 1024^2 },
			"GB" = { denominator <- 1024^3 },
			"TB" = { denominator <- 1024^4 }
	)
	
	ram <- nNodes^2*double.bytes*nCopies / denominator
	
	return(list("RAM" = ram, "units" = units))
}
