print.set <- 
function(x, ...) {
### summarising data set object information
	
	cat("\n\tDataset\n\n")
	if(!is.null(x$description)) cat("description:", x$description, "\n")
	cat("samples:", nrow(x$data), "\n")
	cat("height:", x$height, attr(x$height, "unit"), "\n")
	signals <- names(x$data)
	if(!is.null(x$data$v.avg) && !is.null(x$data$dir.avg)) cat(paste0("availability: ", round(sum(!is.na(x$data$v.avg) & !is.na(x$data$dir.avg)) * 100 / nrow(x$data)), "%\n")) else cat("availability: 0%\n")
	cat("\nsignals:")
	if(is.null(attr(x$data, "clean"))) clean <- "not cleaned"
	else clean <- unlist(attr(x$data, "clean"))
	sig.tbl <- data.frame(matrix(NA, nrow=length(signals), ncol=1))
	names(sig.tbl) <- " "
	row.names(sig.tbl) <- signals
	for(i in 1:length(signals)) if(any(signals==row.names(sig.tbl)[i])) sig.tbl[i,1] <- "(original)"
	if(length(clean)!=1 && clean[1]!="not cleaned") {
		if(any(names(clean)=="v.avg.min" || names(clean)=="v.avg.max")) if(any(signals=="v.avg")) sig.tbl[which(row.names(sig.tbl)=="v.avg"),1] <- "(cleaned) "
		if(any(names(clean)=="dir.clean")) if(clean[names(clean)=="dir.clean"]) if(any(signals=="dir.avg")) sig.tbl[which(row.names(sig.tbl)=="dir.avg"),1] <- "(cleaned) "
		if(any(names(clean)=="icing")) if(clean[names(clean)=="icing"]) if(any(signals=="dir.avg")) sig.tbl[which(row.names(sig.tbl)=="dir.avg"),2] <- "(cleaned) "
		if(any(names(clean)=="turb.clean")) if(any(signals=="turb.int")) sig.tbl[which(row.names(sig.tbl)=="turb.int"),1] <- "(cleaned) "
	}
	sig.tbl[is.na(sig.tbl)] <- ""
	print(sig.tbl, quote=FALSE)
}
