`fdbkFunc` <-
function(y) {
	if(is.list(y) != TRUE) {
		stop("Error: ", paste(sQuote("y"), sep = ""), " should be a list.")
	}
	u <- y
	P <- nrow(y[[1]])
	T <- ncol(y[[1]])
	for(i in 1:length(y)) {
		u[[i]] <- cbind(rep(0, P), u[[i]][,-T])
	}
	return(u)
}

