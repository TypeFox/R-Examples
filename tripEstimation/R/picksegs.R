

## pick out ranges from data - an even number length vector
pick <- function(id, val, nsee = 10000) {

	ini <- 1:nsee
	
	n <- length(id)
	if (n != length(val)) stop("data of different lengths")
	
	nn <- n%/%nsee + 1
	
	cutsub <- cut(1:n, nn)
	
	picks <- NULL
	for (i in 1:nn) {
		sub <- cutsub == levels(cutsub)[i]
		plot(id[sub], val[sub], pch = ".")
		picks <- c(picks, locator()$x	)	
	}
	
	return(picks)
}


## generate segments based on the pairs 
picksegs <- function(twind, n) {
	seg <- 0
	segments <- rep(NA, n) 
	twind <- round(twind)
	for (i in seq(1, length(twind), by = 2)) {
		seg <- seg + 1
		rg <- seq(twind[i], twind[i+1]) 
		segments[rg] <- seg
	}
	segments
}