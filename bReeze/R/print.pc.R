print.pc <- 
function(x, ...) {
### summarising power curve object information
	
	#cat("\n\tPower curve", substitute(x), "\n\n")
	cat("\n\tPower curve\n\n")
	if(!is.null(attr(x, "description"))) cat("description:", attr(x, "description"), "\n")
	cat("rated power:", attr(x, "rated.power"), "\n")
	cat("air pressure:", attr(x, "rho"), "\n\n")
	tbl.units <- data.frame(t(names(x)))
	tbl.units[,] <- "[-]"
	tbl.units[,1] <- paste0("[", attr(x, "units")[1], "]")
	tbl.units[,2] <- paste0("[", attr(x, "units")[2], "]")
	x[is.na(x)] <- ""
	obj <- as.data.frame(lapply(x, as.character))
	names(x)[1] <- "wind speed"
	names(x)[2] <- "power"
	names(x)[names(x)=="cp"] <- "power coefficient"
	names(x)[names(x)=="ct"] <- "thrust coefficient"
	names(tbl.units) <- names(obj) <- names(x)
	row.names(tbl.units) <- " "
	row.names(obj) <- as.character(1:nrow(obj))
	print(rbind(tbl.units, obj), quote=FALSE)
	cat("\n")
}
