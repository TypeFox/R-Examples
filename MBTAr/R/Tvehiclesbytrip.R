Tvehiclesbytrip <- function(trip_id,api_key){
  # finds vehicle locations for a particular trip
  query <- "vehiclesbytrip"
  base_url <- paste("http://realtime.mbta.com/developer/api/v2/",query,"?api_key=",api_key,"&format=json",sep="")
  full_url <- paste(base_url,"&trip=",trip_id,sep="")
  if(length(grep("Green",x=Tschedulebytrip(trip_id=trip_id,api_key=api_key)$route_id))>0){
    stop("Error: Predictions and real time tracking not yet available for the Green Line.")
  }
  rawdata <- readLines(full_url, warn = F) # add error if 404 not found - warn that trip must be in the next hour.
  dl <- jsonlite::fromJSON(txt=rawdata,simplifyDataFrame = T)
  allout <- data.frame(
    route_id=dl$route_id,
    route_name=dl$route_name,
    route_type=dl$route_type,
    mode_name=dl$mode_name,
    trip_id=dl$trip_id,
    trip_name=dl$trip_name,
    trip_headsign=dl$trip_headsign,
    direction_id=dl$direction_id,
    direction_name=dl$direction_name,
    vehicle_id=dl$vehicle$vehicle_id,
    vehicle_lat=dl$vehicle$vehicle_lat,
    vehicle_lon=dl$vehicle$vehicle_lon,
    vehicle_bearing=ifelse(length(dl$vehicle$vehicle_bearing)>0,dl$vehicle$vehicle_bearing,NA) # not available for bus
  )
  return(allout)
}