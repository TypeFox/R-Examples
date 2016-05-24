second.by.second <-
function(data)
{
	sec.by.sec.data <- data.frame(time=NA, date=NA, ap.posture=NA, mets=NA, met.hours=NA, steps=NA)
  sec.by.sec.data <- sec.by.sec.data[-1,]
  
  #	AP assigns 1.4 METs to standing.  We believe standing is a light intensity activity so we change this value to 1.5 METs in order for standing to be counted towards light intensity time and to adjust to MET-hrs.
		data$interval <- as.numeric(data$interval)
		
		data$methrs[data$activity==1] <- (data$interval[data$activity==1]/3600)*1.5
		
		data$methrs <- as.numeric(data$methrs)
  
  	n <- dim(data)[1]
  	time.of.each.event <- as.vector(difftime(strptime(data$time[seq_len(n - 1) + 1],format="%Y-%m-%d %H:%M:%S"),strptime(data$time[seq_len(n - 1)],format="%Y-%m-%d %H:%M:%S"), units="secs"))
	start.time <- strptime(data$time[1],format="%Y-%m-%d %H:%M:%S")

  	time.of.each.event <- c(time.of.each.event, round(data[n,"interval"],0))
	te <- length(time.of.each.event)
	time.of.each.event[is.na(time.of.each.event)==T] <- 1
	events <- rep((1:te),time.of.each.event)
	
	acts <- rep(data$activity,time.of.each.event)
	n <- length(acts)
# The met hours per second in the interval.
	met.hours <- data$methrs/data$interval 	
	met.hours <- rep(met.hours,time.of.each.event)
# To compute mets per second in the interval, multiply methours by 3600 sec/hour and divide by number of seconds.
		mets <- data$methrs * 3600/data$interval
		mets <- rep(mets,time.of.each.event)
		steps <- rep(data$cumulativesteps,time.of.each.event)
# Make 15-sec epoch variable and METs
		times <- start.time+(0:(n-1))
		fifteen.sec.times <- start.time + (15*rep(0:(floor(n/15)),each=15,length=n))
		fifteen.sec.mets <- tapply(mets, fifteen.sec.times, mean)
		fifteen.sec.mets <- rep(fifteen.sec.mets, each=15, length=n)

# Make 1-min epoch variable and METs
		times <- start.time+(0:(n-1))
		one.min.times <- start.time + (60*rep(0:(floor(n/60)),each=60,length=n))
		one.min.mets <- tapply(mets, one.min.times, mean)
		one.min.mets <- rep(one.min.mets, each=60, length=n)
		
		date <- substring(format(times),1,10)

sec.by.sec.data <- merge(sec.by.sec.data, data.frame(time=times, date=date, ap.posture=acts, mets=mets, fifteen.sec.mets=fifteen.sec.mets, one.min.mets=one.min.mets, met.hours=met.hours, steps=steps, num.events=events, stringsAsFactors=FALSE), all=TRUE)

sec.by.sec.data$mets <- signif(sec.by.sec.data$mets,3)

return(sec.by.sec.data)
}

