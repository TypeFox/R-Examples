wearingMarking <-
function(dataset, 
                          frame = 90, 
                          perMinuteCts = 60,
                          TS = "TimeStamp",
                          cts = "counts", 
                          streamFrame = NULL, 
                          allowanceFrame= 2, 
                          newcolname = "wearing",
                          getMinuteMarking = FALSE,
                          dayStart = "00:00:00",
                          dayEnd = "23:59:59",
                          ...)
{
    if(perMinuteCts != 1){
        #not a minute data run collapse
        data2 = dataCollapser(dataset, TS=TS, by = 60, col = cts, ...)
    }else{
        data2 = dataset
    }
    data3 = marking(data2, frame = frame, cts = cts, streamFrame = streamFrame, 
                          allowanceFrame = allowanceFrame, newcolname = newcolname)

    colName = names(data3)
    if(!getMinuteMarking){
        dataset$key = substring(dataset[,names(dataset)[TS ==  names(dataset)]], 1, 16)
        data3$key = substring(data3[,names(data3)[TS ==  names(data3)]], 1, 16)
        data4 = merge(dataset, data3[c(newcolname, "key")], all.x = TRUE, by = "key")[c(colName)]
    }else{
        data4 = data3[c(colName)]
    }

    data4$weekday = weekdays(as.POSIXlt(as.vector(data4[,TS]),format = "%Y-%m-%d %H:%M:%S", tz = "GMT"))
    markingTime(data4, TS, dayStart, dayEnd)
}

