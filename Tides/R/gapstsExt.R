gapsts <- function(ts,            	# array of times, consisting of different continuous subseries seperated by large gaps
                    dtMax,        	# maximum time interval in a continuous series, or equivalently minimum interval to be characterized as gap.
                    unit = "mins"  	# unit of dtMax; used when ts belongs to class POSIXt
                    )
{

if (!inherits(ts,"POSIXt")){
timediffs <- ts[1:(length(ts)-1)] - ts[2:(length(ts))]
} else {
timediffs <- difftime(ts[1:(length(ts)-1)], ts[2:(length(ts))],units=unit)
}

if (!any(timediffs < - dtMax)) return(NULL)
#Select gaps > dtMax in a timeseries ts
gaps <- ts[c(timediffs < -dtMax,FALSE)]
gaps <- data.frame(t1 = gaps)
gaps$t2 <- ts[match(gaps$t1,ts)  + 1]
gaps$n <- 1:dim(gaps)[1]

if (!inherits(ts,"POSIXt")){
gaps$dt <- gaps$t2 - gaps$t1
} else {
gaps$dt <- difftime(gaps$t2,gaps$t1,units=unit)
}

return(gaps) #Data frame with the initial time, end time and time difference (unit = unit) of each interval > dtMax
}



