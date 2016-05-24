`assemble.date` <-
function (table) 
{
    t.midnight <- table$Hour==24
    table$Hour[t.midnight] <- 23 
    dates <- strptime(paste(sep = "", table$Year, "-", table$Month,
        "-", table$Day, " ", table$Hour%%24), "%Y-%m-%d %H") 
    dates[t.midnight] <- dates[t.midnight]+3600
    return(as.POSIXct(dates,tz="GMT"))
}
