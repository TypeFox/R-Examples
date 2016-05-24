modified <-
function (dat, d1, d2) 
{
    da <- as.POSIXct(d1, format = "%d-%m-%Y %H:%M:%S", tz = Sys.timezone())
    db <- as.POSIXct(d2, format = "%d-%m-%Y %H:%M:%S", tz = Sys.timezone())
    d8 <- as.POSIXct(dat$Modified, format = "%d-%m-%Y %H:%M:%S", 
        tz = Sys.timezone())
    f1 <- which(d8 >= da & d8 <= db)
    f2 <- order(d8[f1], decreasing = T)
    f <- f1[f2]
}
