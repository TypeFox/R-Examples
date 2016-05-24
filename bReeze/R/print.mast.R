print.mast <- 
function(x, ...) {
### summarising met mast object information
	
	if(!is.null(x$location)) loc <- x$location
	if(!is.null(x$description)) desc <- x$description
	num.sets <- length(x$sets)
	num.samples <- nrow(x$sets[[1]]$data)
	heights <- x$sets[[1]]$height
	h.unit <- attr(x$sets[[1]]$height, "unit")
	interval <- x$timestamp[2]-x$timestamp[1]
	if(attr(interval, "units")=="days" && interval>1) warning("Availability cannot be calculated - time interval longer than 1 day", call.=FALSE)
	if(attr(interval, "units")=="days") daily.samples <- interval
	if(attr(interval, "units")=="hours") daily.samples <- 24/as.numeric(interval)
	if(attr(interval, "units")=="mins") daily.samples <- 24*60/as.numeric(interval)
	if(attr(interval, "units")=="secs") daily.samples <- 24*60*60/as.numeric(interval)
	period.start <- as.character(x$timestamp[1])
	period.end <- as.character(x$timestamp[num.samples])
	if(nchar(period.start)==10) period.start <- paste(period.start, "00:00:00")
	if(nchar(period.end)==10) period.end <- paste(period.end, "00:00:00")
	period.start <- strptime(period.start, format="%Y-%m-%d %H:%M:%S")
	period.end <- strptime(period.end, format="%Y-%m-%d %H:%M:%S")
	period.days <- as.numeric(period.end-period.start)
	signals <- names(x$sets[[1]]$data)
	if(is.null(x$sets[[1]]$data$v.avg)) wind.speed <- 0
	else wind.speed <- mean(x$sets[[1]]$data$v.avg, na.rm=TRUE); v.unit <- attr(x$sets[[1]]$data$v.avg, "unit")
	if(is.null(x$sets[[1]]$data$v.avg) || is.null(x$sets[[1]]$data$dir.avg)) avail <- 0 
	else avail <- sum(!is.na(x$sets[[1]]$data$v.avg) & !is.na(x$sets[[1]]$data$dir.avg)) * 100 / (daily.samples*period.days)
	if(is.null(attr(x$sets[[1]]$data, "clean"))) clean <- "not cleaned"
	else clean <- unlist(attr(x$sets[[1]]$data, "clean"))
	
	if(num.sets>1) {
		signals <- list(signals)
		clean <- list(clean)
		for(i in 2:num.sets) {
			heights <- append(heights, x$sets[[i]]$height)
			signals[[i]] <- names(x$sets[[i]]$data)
			if(is.null(x$sets[[i]]$data$v.avg)) wind.speed <- 0 
			else wind.speed <- append(wind.speed, mean(x$sets[[i]]$data$v.avg, na.rm=TRUE)); v.unit <- attr(x$sets[[i]]$data$v.avg, "unit")
			if(is.null(x$sets[[i]]$data$v.avg) || is.null(x$sets[[i]]$data$dir.avg)) avail <- append(avail, 0)
			else avail <- append(avail, sum(!is.na(x$sets[[i]]$data$v.avg) & !is.na(x$sets[[i]]$data$dir.avg)) * 100 / (daily.samples*period.days))
			if(is.null(attr(x$sets[[i]]$data, "clean"))) clean[[i]] <- "not cleaned"
			else clean[[i]] <- unlist(attr(x$sets[[i]]$data, "clean"))
		}
	}
	
	attr(heights, "unit") <- h.unit
	attr(avail, "unit") <- "%"
	attr(wind.speed, "unit") <- v.unit
	
	cat(paste("\n\tMet mast\n\n"))
	if(!is.null(x$location)) {
		if(x$location[1]<0) ns <- " degree South, " else ns <- " degree North, "
		if(x$location[2]<0) we <- " degree West" else we <- " degree East"
		cat("location: ", abs(x$location[1]), ns, abs(x$location[2]), we, "\n", sep="")
	}
	if(!is.null(x$description)) cat("description:", x$description, "\n\n")
	cat(paste0("measuring period: from ", period.start, " to ", period.end, " (", round(period.days, 1), " days)\n"))
	cat("samples:", num.samples, "\n\n")
	cat("datasets (", num.sets, "):\n", sep="")
	det <- data.frame(cbind(heights, round(wind.speed, 2), round(avail, 1)))
	tbl.units <- data.frame(t(names(det)))
	tbl.units[,1] <- paste0("[", h.unit, "]")
	tbl.units[,2] <- paste0("[", v.unit, "]")
	tbl.units[,3] <- "[%]"
	det[is.na(det)] <- ""
	det <- as.data.frame(lapply(det, as.character))
	names(det) <- names(tbl.units) <- c("height", "wind speed", "availability")
	row.names(tbl.units) <- " "
	row.names(det) <- names(x$sets)
	print(rbind(tbl.units, det), quote=FALSE)
	sig <- unique(unlist(signals))
	cat("\nsignals (", length(sig), "):\n", sep="")
	sig.tbl <- data.frame(matrix(NA, nrow=length(sig), ncol=num.sets))
	row.names(sig.tbl) <- sig
	names(sig.tbl) <- c(names(x$sets))
	if(num.sets==1) {
		signals <- list(signals)
		clean <- list(clean)
	}
	for(i in 1:length(sig)) for(j in 1:num.sets) if(any(signals[[j]]==row.names(sig.tbl)[i])) sig.tbl[i,j] <- "o"
	for(j in 1:num.sets) {
		if(length(clean[[j]])!=1 && clean[[j]][1]!="not cleaned") {
			if(any(names(clean[[j]])=="v.avg.min" || names(clean[[j]])=="v.avg.max")) if(any(signals[[j]]=="v.avg")) sig.tbl[which(row.names(sig.tbl)=="v.avg"),j] <- "c"
			if(any(names(clean[[j]])=="dir.clean")) if(clean[[j]][names(clean[[j]])=="dir.clean"]) if(any(signals[[j]]=="dir.avg")) sig.tbl[which(row.names(sig.tbl)=="dir.avg"),j] <- "c"
			if(any(names(clean[[j]])=="icing")) if(clean[[j]][names(clean[[j]])=="icing"]) if(any(signals[[j]]=="dir.avg")) sig.tbl[which(row.names(sig.tbl)=="dir.avg"),j] <- "c"
			if(any(names(clean[[j]])=="turb.clean")) if(any(signals[[j]]=="turb.int")) sig.tbl[which(row.names(sig.tbl)=="turb.int"),j] <- "c"
		}
	}
	sig.tbl[is.na(sig.tbl)] <- ""
	print(sig.tbl, quote=FALSE)
	cat("\t(o=original data, c=cleaned data)\n")
}
