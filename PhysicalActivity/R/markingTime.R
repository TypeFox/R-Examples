`markingTime` <-
function(dataset, timestamp, startTime = "00:00:00", endTime = "23:59:59")
{
    if(is.numeric(timestamp)){
        cadval = as.vector(dataset[,timestamp])
    }else {
        cadval = as.vector(dataset[,c(names(dataset)== timestamp)])
    }
    size = length(cadval)

    daystart = paste(substring(as.POSIXlt(cadval[1], tz = "GMT"),1,10), startTime)
    dayend = paste(substring(as.POSIXlt(cadval[1], tz = "GMT"),1,10), endTime)

    days = 1;
    dayMarking = rep(NA, size)
    while(as.POSIXlt(dayend, tz = "GMT")  < 60*60*24 + as.POSIXlt(cadval[size], tz = "GMT"))
    {
        dayMarking[as.POSIXlt(cadval, tz = "GMT")>= as.POSIXlt(daystart, tz = "GMT") & 
                   as.POSIXlt(cadval, tz = "GMT")<= as.POSIXlt(dayend, tz = "GMT")] = days
        
        days = days+1
        daystart = as.POSIXlt(dayend, tz = "GMT")+ 1 
        dayend = as.POSIXlt(dayend, tz = "GMT")+ 60*60*24
    }
    temp = cbind(dataset, days = dayMarking) 
    return(temp)
}

