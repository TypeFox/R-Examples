#' This function retrieves time series of gauged daily flow from the NRFA database.
#'
#' @author Claudia Vitolo
#'
#' @description This function retrieves time series of gauged daily flow from the NRFA database and stores it in zoo format.
#'
#' @param myList this is a nested list which comes from the conversion of waterml2 data.
#'
#' @return time series of class zoo
#'

FindTS <- function(myList){

  #require(zoo)
  #require(stringr)

  temp <- myList[[1]][[1]]$Collection$observationMember$OM_Observation$result$MeasurementTimeseries

  start <- 3
  end <- length(temp) - 1

  temp2 <- temp[start:end]
  l <- length(temp2)

  time <- c()
  value <- c()

  for (x in 1:l) {

    time[x] <- str_sub(temp[x]$point$MeasurementTVP$time,start=2,end=-1)
    value[x] <- as.numeric(temp[x]$point$MeasurementTVP$value)

  }

  myTS <- zoo(value, order.by=as.POSIXlt(time))

  if ( any(is.na(index(myTS))) ){

    toDelete <- which(is.na(index(myTS)))
    myTS <- myTS[-toDelete]

  }

  return(myTS)

}
