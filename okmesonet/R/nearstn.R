nearstn <- function(pnt.lon, pnt.lat, startdate, enddate) {
  ## Calculate the nearest OK Mesonet station given latitude and longitude
  ## coordinates.
  ## Retrieves station coordinates from stations data frame
  ## 
  ## Arguments:
  ##  pnt.lon: longitude of location, given in decimal degrees
  ##  pnt.lat: latitude of location, given in decimal degrees
  ##  startdate:  start date of desired data
  ##  enddate: end date of desired data
  ## Returns: four letter station identifier as character object
  
  ## subset active stations for given time frame
  stations.active <- subset(okstations, Commissioned<=startdate & 
                            Decommissioned>enddate)
  
  ## calculate station distances from point location
  stndistance <- mapply(FUN=vincenty, stations.active$Longitude,
                        stations.active$Latitude,
                        MoreArgs=list(lon1=pnt.lon, lat1=pnt.lat))
  
  ## determine nearest station
  nearstn <- which.min(stndistance)
  ## distance in meters
  neardist <- stndistance[nearstn]
  
  ## print message displaying station
  msg <- paste("Using ", stations.active$Name[nearstn], " (", 
               stations.active$Identifier[nearstn], ") station. ", 
               "Distance: ", neardist/1000, " km ",
               "Commissioned: ", stations.active$Commissioned[nearstn], 
               " Decommissioned: ", stations.active$Decommissioned[nearstn], 
               sep="")
  message(msg)
  
  ## return four letter station identifier as lowercase
  return(tolower(stations.active$Identifier[nearstn]))
}