modifiedtoday <-
function (dat) 
{
    d8 <- as.POSIXct(dat$Modified, format = "%d-%m-%Y %H:%M:%S", 
        tz = Sys.timezone())
    tod <- format(Sys.time(), "%d-%m-%Y")
    tod <- paste(tod, "00:00:00", sep = " ")
    today <- as.POSIXct(tod, format = "%d-%m-%Y %H:%M:%S", tz = Sys.timezone())
    f1 <- which(d8 > today)
    f2 <- order(d8[f1], decreasing = T)
    f <- f1[f2]
}
