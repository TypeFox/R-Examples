`dataCollapser` <-
function(dataset, TS, by, col, func = sum, ...)
{
    ts = as.vector(dataset[,TS])
    ct = as.numeric(dataset[,col])

    timeRange = range(as.vector(ts))
    epoch = as.numeric(as.POSIXlt(ts[2]) - as.POSIXlt(ts[1]))
    ratio = by/epoch

    newrange = c(0: (ceiling(length(ts)/ratio)-1))*by
    step1 = rep(as.POSIXlt(timeRange[1], tz = "GMT"), length(newrange))
    newts = gsub(" GMT", "", step1+ newrange)
    newct = rep(NA, length(newrange))

    i = 1
    prevts = 
    while(i <= length(newts))
    {
        start = (i-1)*ratio +1
        end = i*ratio
        if(end > length(ct))
        {    end = length(ct)} 

        newct[i] = func(ct[start:end], ...)
        
        i = i+1
    }

    tf = data.frame(timestamp = newts, counts = newct)
    names(tf) = c(TS, col)
    return(tf)
}

