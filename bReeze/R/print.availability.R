print.availability <- 
function(x, ...) {
### summarising availability object information
	
	cat("\n\tAvailability for pairs of wind speed and direction\n\n")
	tot <- x[[1]]$total
	if(length(x)>1) for(i in 2:length(x)) tot <- rbind(tot, x[[i]]$total)
	tbl.units <- data.frame(t(names(tot)))
	tbl.units[,] <- "[d]"
	tbl.units[,1] <- "[%]"
	names(tot)[2:3] <- c("effective period", "total period")
	names(tbl.units) <- names(tot)
	row.names(tbl.units) <- " "
	row.names(tot) <- names(x)
	print(rbind(tbl.units, tot), quote=FALSE)
	cat("\nnumber of daily samples:\n")
	cat(names(x)[1], "\n")
	x[[1]]$daily[is.na(x[[1]]$daily)] <- ""
	names(x[[1]]$daily)[1] <- "%"
	print(x[[1]]$daily, quote=FALSE)
	cat("\n")
	if(length(x)>1) {
		for(i in 2:length(x)) {
			cat(names(x)[i], "\n")
			x[[i]]$daily[is.na(x[[i]]$daily)] <- ""
			names(x[[i]]$daily)[1] <- "%"
			print(x[[i]]$daily, quote=FALSE)
			cat("\n")
		}
	}
	if(attr(x, "call")$v.set[1]=="all") attr(x, "call")$v.set[1] <- "\"all\""
	if(attr(x, "call")$dir.set[1]=="all") attr(x, "call")$dir.set[1] <- "\"all\""
	if(length(attr(x, "call")$v.set)==1) vset <- paste0(", v.set=", attr(x, "call")$v.set)
	else vset <- paste0(", v.set=c(", paste0(attr(x, "call")$v.set, collapse=", "), ")")
	if(length(attr(x, "call")$dir.set)==1) dirset <- paste0(", dir.set=", attr(x, "call")$dir.set)
	else dirset <- paste0(", dir.set=c(", paste0(attr(x, "call")$dir.set, collapse=", "), ")")
	if(!any(!is.na(attr(x, "call")$subset))) subs <- ", subset=NA"
	else subs <- paste0(", subset=c(\"", paste(attr(x, "call")$subset, collapse="\", \""), "\")")
	cat("call: availability(mast=", attr(x, "call")$mast, vset, dirset, subs, ", digits=", attr(x, "call")$digits, ", print=", attr(x, "call")$print, ")\n\n", sep="")
}
