print.month.stats <- 
function(x, ...) {
### summarising month.stats object information
	
	cat("\n\tMonthly statistics\n\n")
	cat(names(x)[1], "\n")
	x[[1]][is.na(x[[1]])] <- ""
	if(length(row.names(x[[1]]))==14) names(x[[1]])[length(names(x[[1]]))] <- row.names(x[[1]])[14] <- gsub("\\.", " ", row.names(x[[1]])[14])
	row.names(x[[1]])[1:12] <- c(toupper(row.names(x[[1]])[1:12]))
	print(x[[1]], quote=FALSE)
	cat("\n")
	if(length(x)>1) {
		for(i in 2:length(x)) {
			cat(names(x)[i], "\n")
			x[[i]][is.na(x[[i]])] <- ""
			if(length(row.names(x[[1]]))==14) names(x[[i]])[length(names(x[[i]]))] <- row.names(x[[i]])[14] <- gsub("\\.", " ", row.names(x[[i]])[14])
			row.names(x[[i]])[1:12] <- c(toupper(row.names(x[[i]])[1:12]))
			print(x[[i]], quote=FALSE)
			cat("\n")
		}
	}
	if(attr(x, "call")$set=="all") attr(x, "call")$set <- "\"all\""
	if(!any(!is.na(attr(x, "call")$subset))) subs <- ", subset=NA"
	else subs <- paste0(", subset=c(\"", paste(attr(x, "call")$subset, collapse="\", \""), "\")")
	cat("call: month.stats(mast=", attr(x, "call")$mast, ", set=", attr(x, "call")$set, ", signal=\"", attr(x, "call")$signal, "\", fun=\"", attr(x, "call")$fun, "\"", subs, ", digits=", attr(x, "call")$digits, ", print=", attr(x, "call")$print, ")\n\n", sep="")
}
