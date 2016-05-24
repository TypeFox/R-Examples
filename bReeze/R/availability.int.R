availability.int <-
function(v.avg, dir.avg, ts, start.year, start.month, num.months, period.days, digits) {
### internal function for calculation of availability

	num.samples <- length(ts)
	avail <- data.frame(matrix(0, ncol=32, nrow=num.months))
	names(avail) <- c("all", 1:31)
	yr <- start.year
	mon <- start.month
	for(i in 1:num.months) {
		if(mon<10) row.names(avail)[i] <- paste(yr, paste0("0", mon), sep="-")
		if(mon>=10) row.names(avail)[i] <- paste(yr, mon, sep="-")
		mon <- mon+1
		if(mon==13) {
			yr <- yr+1
			mon <- 1
		}
	}
	
	interval <- ts[2]-ts[1]
	if(attr(interval, "units")=="days" && interval>1) stop("Time interval longer than 1 day")
	if(attr(interval, "units")=="days") daily.samples <- interval
	if(attr(interval, "units")=="hours") daily.samples <- 24/as.numeric(interval)
	if(attr(interval, "units")=="mins") daily.samples <- 24*60/as.numeric(interval)
	if(attr(interval, "units")=="secs") daily.samples <- 24*60*60/as.numeric(interval)
	
	for(m in 1:dim(avail)[1]) {
		yr <- as.numeric(strsplit(row.names(avail)[m], "-")[[1]][1])
		mon <- as.numeric(strsplit(row.names(avail)[m], "-")[[1]][2])
		for(d in 2:32) {
			avail[m,d] <- length(v.avg[ts$year==yr-1900 & ts$mon==mon-1 & ts$mday==d-1 & !is.na(v.avg) & !is.na(dir.avg)])
		}
		if(any(mon==c(1,3,5,7,8,10,12))) days <- 31
		if(any(mon==c(4,6,9,11))) days <- 30
		leap.year <- FALSE
		if(yr%%4==0) {
			leap.year <- TRUE
			if(yr%%100==0) {
				leap.year <- FALSE
				if(yr%%400==0) leap.year <- TRUE
			}
		}
		if(mon==2 & leap.year) days <- 29
		if(mon==2 & !leap.year) days <- 28
		if(days<31) avail[m,(days+2):32] <- NA
		avail[m,1] <- round(sum(avail[m,2:(days+1)]) * 100 / (days*daily.samples), digits)
	}
	
	availability <- sum(!is.na(v.avg) & !is.na(dir.avg)) * 100 / (period.days*daily.samples)
	total <- data.frame(availability=round(availability, digits), effective.period=round(availability / 100 * period.days, digits), total.period=round(period.days, digits))
	
	attr(avail, "num.daily.samples") <- daily.samples
	attr(total$availability, "unit") <- "%"
	attr(total$effective.period, "unit") <- "d"
	attr(total$total.period, "unit") <- "d"
	
	return(list(total=total, daily=avail))
}
