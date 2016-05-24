print.energy <- 
function(x, ...) {
### summarising energy object information
	
	cat("\n\tWind energy content\n\n")
	obj <- x[[1]]
	if(length(x)>1) {
		for(i in 2:length(x)) {
			obj <- rbind(obj, x[[i]])
		}
		obj <- as.data.frame(t(as.matrix(obj)))
	}	
	obj <- as.data.frame(obj)
	names(obj) <- names(x)
	row.names(obj) <- c(toupper(head(attr(x, "row.names"), -1)), tail(attr(x, "row.names"), 1))
	obj[obj==0] <- ""
	print(obj, quote=FALSE)
	cat("\t(all values in ", attr(x, "unit"), ")\n", sep="")
	cat("\ncall: energy(wb=", attr(x, "call")$wb, ", rho=", attr(x, "call")$rho, ", bins=c(", paste(attr(x, "call")$bins, collapse=", "), "), digits=", attr(x, "call")$digits, ", print=", attr(x, "call")$print, ")\n\n", sep="")
}
