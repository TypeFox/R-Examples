wclim <-
function (f = 1:19, full = F) 
{
    bio <- {
    }
    for (i in 1:9) {
        bio[i] <- paste("BIO0", i, sep = "")
    }
    for (i in 10:19) {
        bio[i] <- paste("BIO", i, sep = "")
    }
    abbrev <- c("AMT", "MDR", "I", "TS", "MTWM", "MTCM", "TAR", 
        "MTWeQ", "MTDQ", "MTWQ", "MTCQ", "AP", "PWeM", "PDM", 
        "PS", "PWeQ", "PDQ", "PWQ", "PCQ")
    name <- c("Annual Mean Temperature", "Mean Diurnal Range", 
        "Isothermality", "Temperature seasonality", "Max Temperature of Warmest Month", 
        "Min Temperature of Coldest Month", "Temperature Annual Range", 
        "Mean Temperature of Wettest Quarter", "Mean Temperature of Driest Quarter", 
        "Mean Temperature of Warmest Quarter", "Mean Temperature of Coldest Quarter", 
        "Annual Precipitation", "Precipitation of Wettest Month", 
        "Precipitation of Driest Month", "Precipitation Seasonality", 
        "Precipitation of Wettest Quarter", "Precipitation of Driest Quarter", 
        "Precipitation of Warmest Quarter", "Precipitation of Coldest Quarter")
    correction <- c(rep(10, 11), rep(1, 8))
    wclim <- data.frame(bio, abbrev, name, correction, stringsAsFactors = F)
    if (full == T) {
        wcm <- wclim[f, ]
    }
    else {
        wcm <- wclim$name[f]
    }
    return(wcm)
}
