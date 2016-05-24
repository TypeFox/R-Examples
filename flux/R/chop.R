## prepare function, making a list with a data.frame for
## each flux measurement (one chamber placement)
chop <-
function(dat, factors, nmes = NULL, min.cm = 3){
	## make the list for the by function
	sellist <- lapply(c(1:length(factors)), function(x) dat[,factors[x]])
	## separate into tables per chamber measurement
	conz.parts <- split(dat, sellist, drop=TRUE)
	## tables with less then a specified number of concentration
	## measurements (via min.cm) are eliminated
	flux.sel <- sapply(conz.parts, function(x) nrow(x)>=min.cm)
	conz.parts <- conz.parts[flux.sel]
	## do the naming (for easy access to the data)
	if(is.null(nmes)){ nmes <- factors }
	nams <- data.frame(t(sapply(conz.parts, function(x) sapply(x[][1,nmes], as.character))))
	nams$all <- apply(nams, 1, function(x) paste(sapply(x, as.character), sep=".", collapse="."))
	names(conz.parts) <- nams$all
	conz.parts <- list(tables = conz.parts, nmes = nams)
	return(conz.parts)
} 