print.profile <- 
function(x, ...) {
### summarising profile object information
	
	cat("\n\tWind profile\n\n")
	tbl.units <- data.frame(t(names(x$profile)))
	tbl.units[,1] <- paste0("[", attr(x$profile, "unit")[1], "]")
	tbl.units[,2] <- paste0("[", attr(x$profile, "unit")[2], "]")
	obj <- as.data.frame(lapply(x$profile, as.character))
	names(x$profile)[2] <- "wind speed"
	names(tbl.units) <- names(obj) <- names(x$profile)
	row.names(tbl.units) <- " "
	row.names(obj) <- c(toupper(head(row.names(x$profile), -1)), tail(row.names(x$profile), 1))
	print(rbind(tbl.units, obj), quote=FALSE)
	cat("\nreference height:", x$h.ref, attr(x$h.ref, "unit"), "\n")
	if(is.null(attr(x, "call")$alpha)) alph <- ", alpha=NULL"
	else {
		if(length(attr(x, "call")$alpha)==1) alph <- paste0(", alpha=", attr(x, "call")$alpha)
		else  alph <- paste0(", alpha=c(", paste0(attr(x, "call")$alpha, collapse=", "), ")")
	}
	if(length(attr(x, "call")$v.set)==1) vset <- paste0(", v.set=", attr(x, "call")$v.set)
	else vset <- paste0(", v.set=c(", paste0(attr(x, "call")$v.set, collapse=", "), ")")
	if(!any(!is.na(attr(x, "call")$subset))) subs <- ", subset=NA"
	else subs <- paste0(", subset=c(\"", paste(attr(x, "call")$subset, collapse="\", \""), "\")")
	cat("\ncall: profile(mast=", attr(x, "call")$mast, vset, ", dir.set=", attr(x, "call")$dir.set, ", num.sectors=", attr(x, "call")$num.sectors, ", method=\"", attr(x, "call")$method, "\"", alph, subs, ", digits=", attr(x, "call")$digits, ", print=", attr(x, "call")$print, ")\n\n", sep="")
}
