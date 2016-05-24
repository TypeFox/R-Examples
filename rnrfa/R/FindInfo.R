#' Extract info from the waterml2.
#'
#' @author Claudia Vitolo
#'
#' @description This function retrieves the most important info from the waterml2 file given as a list.
#'
#' @param myList this is a nested list which comes from the conversion of waterml2 data
#'
#' @return named vector containing the following information: stationName, Latitude, Longitude, typeOfMeasurement, timeZone, remarks
#'

FindInfo <- function(myList){

  #require(stringr)

  x <- myList[[1]][[1]]$Collection

  stationName <- x$observationMember$OM_Observation$featureOfInterest[[2]]

  typeOfMeasurement <- x$localDictionary$Dictionary$dictionaryEntry$Definition$remarks[[1]]
  variable     <- unlist(strsplit(typeOfMeasurement, ","))[[1]]
  units        <- unlist(strsplit(typeOfMeasurement, ","))[[2]]
  typeFunction <- unlist(strsplit(typeOfMeasurement, ","))[[3]]
  timeStep     <- unlist(strsplit(typeOfMeasurement, ","))[[4]]

  temp <- x$samplingFeatureMember$MonitoringPoint$shape$Point$pos$text
  Latitude <- as.numeric(strsplit(temp, " ")[[1]][[1]])
  Longitude <- as.numeric(strsplit(temp, " ")[[1]][[2]])

  timeZone <- x$samplingFeatureMember$MonitoringPoint$timeZone

  remarks <- x$localDictionary$Dictionary[3]$dictionaryEntry$Definition$remarks

  info <- data.frame(cbind(stationName, Latitude, Longitude, variable, units,
                           typeFunction, timeStep, timeZone, remarks))

  return(info)

}
