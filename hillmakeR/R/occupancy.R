# Get data.frame that describes units used by package
getdfAllowedUnits <- function(){
  unitname <- c("sec", "min", "hour", "day")
  conversion <- c(1, 60, 3600, 86400)
  returnFrame <- data.frame(unitname, conversion)
  return(returnFrame)
}

occupancy <- function(startTimes, stopTimes, resolution = "min", initial = NULL, fillup = NULL, countlast = TRUE){
  #library.dynam("updateArray")
  #orginTime<-as.POSIXct('01/01/1970 00:00',format ="%m/%d/%Y %H:%M", tz = "GMT")
  
 # Check in input is in correct format
  if(!identical(class(startTimes), c("POSIXct","POSIXt")))
    stop(paste("startTimes is class ", class(startTimes),". Must be class POSIXct POSIXt.", sep = ""))
  
# Check in output is in correct format
  if(!identical(class(stopTimes), c("POSIXct","POSIXt")))
    stop(paste("startTimes is class ", class(stopTimes),". Must be class POSIXct POSIXt.", sep = ""))

# Check if both initial and fillup are specified
if(!is.null(fillup) & !is.null(initial))
  stop("It does not make sense to det both initial and fillup parameters.")

  dfAllowedUnits <- getdfAllowedUnits()
  # Check if resolution is an allowed unit
    if (!(resolution %in% dfAllowedUnits$unitname)) {
      allowedunitspaste <-  paste(dfAllowedUnits$unitname, collapse = ", ")
    stop(paste('Resolution can only be set to the following:', allowedunitspaste))
  }
  
  # Calculate LOS (length-of-stay) for each item, remove any negatives. units requires an S added to resolution
  LOS <- as.numeric(stopTimes - startTimes, units = paste(resolution, "s", sep = ""))
  badrecords <- which(LOS <= 0 | is.na(LOS))
  if(length(badrecords) > 0) {
    warning(paste(length(badrecords), " removed because duration is either negative or NA", sep = ""))
    startTimes <- startTimes[-badrecords]
    stopTimes<- stopTimes[-badrecords]
  }
  
  earliestTime <-  min(startTimes)
  latestTime <-  max(stopTimes)
  
  # Get conversion factor
  conversion <- dfAllowedUnits[dfAllowedUnits$unitname == resolution, "conversion"]
  
  startMin <- ceiling((as.numeric(startTimes) -  as.numeric(min(startTimes)))/conversion)
  stopMin  <- ceiling((as.numeric(stopTimes) -   as.numeric(min(startTimes)))/conversion)
  totalMin <- length(startMin)
  
  # consider if the last time unit should be counted (i.e., count the minute that something departs?)
  if(countlast){
    countlastpassed <- 0
  } else {
    countlastpassed <- 1
  }

  counts<-.C("updateArray", size = as.integer(totalMin), starts = as.double(startMin), stops = as.double(stopMin), returnX = as.double(rep(0,max(stopMin) + 1)), countlast = as.integer(countlastpassed))$returnX
  
  times <- seq(from = trunc(earliestTime, units = resolution), to = trunc(latestTime,units = resolution), by = resolution)
  returnFrame <- as.data.frame(times)
  returnFrame <- cbind(returnFrame, counts)
  
  #add initial count if specified
  if(!is.null(initial))
    returnFrame$counts <- returnFrame$counts + initial

  #returnFrame$hour <- as.POSIXlt(returnFrame$dateHourList)$hour
  if(!(is.null(fillup))){
    # recalculate LOS (w/o removed records)
    LOS <- as.numeric(stopTimes - startTimes, units = paste(resolution, "s", sep = ""))
    timetoremove <- quantile(LOS, fillup)
    timetoremove <- ceiling(timetoremove) # round if there is a decimal
    warning(paste(timetoremove, " time steps removed based on fillup quantile"))
    returnFrame <- returnFrame[(timetoremove + 1):nrow(returnFrame),]
  }
  
  return(returnFrame)
}