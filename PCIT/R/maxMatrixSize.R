maxMatrixSize <- function(ram, units=c("MB", "bytes", "KB", "GB", "TB"), nCopies=3) {
	units <- match.arg(units)
	switch(units,
			"bytes" = { units <- 0 },
			"KB" = {units <- 1 },
			"MB" = { units <- 2 },
			"GB" = { units <- 3 },
			"TB" = { units <- 4 }
	)
	
	return(floor(sqrt(ram*(1024^units)/(8*nCopies))))
	
}
