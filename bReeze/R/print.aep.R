print.aep <- 
function(x, ...) {
### summarising aep object information
	
	cat("\n\tAnnual energy production\n\n")
	tbl.units <- data.frame(t(names(x$aep)))
	tbl.units[,] <- paste0("[", attr(x$aep[,3], "unit"), "]")
	tbl.units[,1] <- paste0("[", attr(x$aep[,1], "unit"), "]")
	tbl.units[,2] <- paste0("[", attr(x$aep[,2], "unit"), "]")
	x$aep[x$aep==0] <- ""
	obj <- as.data.frame(lapply(x$aep, as.character))
	names(x$aep)[1] <- "wind speed"
	names(tbl.units) <- names(obj) <- names(x$aep)
	row.names(tbl.units) <- " "
	row.names(obj) <- c(toupper(head(row.names(x$aep), -1)), tail(row.names(x$aep), 1))
	print(rbind(tbl.units, obj), quote=FALSE)
	cat("\ncapacity factor:", x$capacity, "\n")
	cat("\ncall: aep(profile=", attr(x, "call")$profile, ", pc=", attr(x, "call")$pc, ", hub.h=", attr(x, "call")$hub.h, ", rho=", attr(x, "call")$rho, ", avail=", attr(x, "call")$avail, ", bins=c(", paste(attr(x, "call")$bins, collapse=", "), "), sectoral=", attr(x, "call")$sectoral, ", digits=c(", paste(attr(x, "call")$digits, collapse=", "), "), print=", attr(x, "call")$print, ")\n\n", sep="")
}
