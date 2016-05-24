print.turbulence <- 
function(x, ...) {
### summarising turbulence object information
	
	cat("\n\tTurbulence intensity\n\n")
	obj <- x[[1]]
	if(length(x)>1) for(i in 2:length(x)) obj <- cbind(obj, x[[i]])
	obj <- as.data.frame(obj)
	names(obj) <- names(x)
	row.names(obj) <- c(toupper(head(attr(x, "row.names"), -1)), tail(attr(x, "row.names"), 1))
	print(obj, quote=FALSE)
	if(!any(!is.na(attr(x, "call")$subset))) subs <- ", subset=NA"
	else subs <- paste0(", subset=c(\"", paste(attr(x, "call")$subset, collapse="\", \""), "\")")
	cat("\ncall: turbulence(mast=", attr(x, "call")$mast, ", turb.set=", attr(x, "call")$turb.set, ", dir.set=", attr(x, "call")$dir.set, ", num.sectors=", attr(x, "call")$num.sectors, ", bins=c(", paste(attr(x, "call")$bins, collapse=", "), ")", subs, ", digits=", attr(x, "call")$digits, ", print=", attr(x, "call")$print, ")\n\n", sep="")
}
