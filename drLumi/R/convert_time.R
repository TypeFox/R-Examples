convert_time <- function(x){
    aux <- strsplit(x," ")
    info <- do.call("rbind", aux)
    date <- info[,1]
    hour <- info[,2]
    if(nchar(hour)==5){
        hour <- paste(hour,":00",sep="")
    }
    date <- chron(dates.= date, format=c(dates="d/m/y"))
    hour <- chron(times.= hour, format=c(times="h:m:s"))
    newdate <-  chron(dates. = date, times. = hour, 
                format=c(dates="d/m/y", times="h:m:s"), 
                out.format=c(dates="d/m/y", times="h:m:s")) 
    return(newdate)
}

