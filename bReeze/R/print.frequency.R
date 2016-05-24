print.frequency <- 
function(x, ...) {
### summarising frequency object information
	
	cat("\n\tFrequency\n\n")
	freq <- x[[1]]
	for(i in 2:length(x)) freq <- cbind(freq, x[[i]])
	freq <- as.data.frame(freq)
	row.names(freq) <- attr(x, "row.names")
	names(freq) <- names(x)
	tbl.units <- data.frame(t(names(x)))
	tbl.units[,] <- paste0("[", attr(x, "unit")[2], "]")
	tbl.units[,1] <- paste0("[", attr(x, "unit")[1], "]")
	freq[freq==0] <- ""
	obj <- as.data.frame(lapply(freq, as.character))
	names(freq)[1] <- "wind speed"
	names(tbl.units) <- names(obj) <- names(freq)
	row.names(tbl.units) <- " "
	row.names(obj) <- c(toupper(head(attr(x, "row.names"), -1)), tail(attr(x, "row.names"), 1))
	print(rbind(tbl.units, obj), quote=FALSE)
	if(!any(!is.na(attr(x, "call")$subset))) subs <- ", subset=NA"
	else subs <- paste0(", subset=c(\"", paste(attr(x, "call")$subset, collapse="\", \""), "\")")
	cat("\ncall: frequency(mast=", attr(x, "call")$mast, ", v.set=", attr(x, "call")$v.set, ", dir.set=", attr(x, "call")$dir.set, ", num.sectors=", attr(x, "call")$num.sectors, ", bins=c(", paste(attr(x, "call")$bins, collapse=", "), ")", subs, ", digits=", attr(x, "call")$digits, ", print=", attr(x, "call")$print, ")\n\n", sep="")
}
