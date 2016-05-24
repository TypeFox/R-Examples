print.weibull <- 
function(x, ...) {
### summarising weibull object information
		
	cat("\n\tWeibull parameters\n\n")
	tbl.units <- data.frame(t(names(x)))
	tbl.units[,1] <- paste0("[", attr(x, "unit")[1], "]")
	tbl.units[,2] <- paste0("[", attr(x, "unit")[2], "]")
	tbl.units[,3] <- paste0("[", attr(x, "unit")[3], "]")
	tbl.units[,4] <- paste0("[", attr(x, "unit")[4], "]")
	obj <- as.data.frame(lapply(x, as.character))
	names(x)[3] <- "wind speed"
	names(tbl.units) <- names(obj) <- names(x)
	row.names(tbl.units) <- " "
	row.names(obj) <- c(toupper(head(attr(x, "row.names"), -1)), tail(attr(x, "row.names"), 1))
	print(rbind(tbl.units, obj), quote=FALSE)
	if(!any(!is.na(attr(x, "call")$subset))) subs <- ", subset=NA"
	else subs <- paste0(", subset=c(\"", paste(attr(x, "call")$subset, collapse="\", \""), "\")")
	cat("\ncall: weibull(mast=", attr(x, "call")$mast, ", v.set=", attr(x, "call")$v.set, ", dir.set=", attr(x, "call")$dir.set, ", num.sectors=", attr(x, "call")$num.sectors, subs, ", digits=", attr(x, "call")$digits, ", print=", attr(x, "call")$print, ")\n\n", sep="")
}
