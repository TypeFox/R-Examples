Tvehiclesbyroute <- function(route_id,api_key){
  # vehicle positions for ongoing and upcoming trips for a particular route
  query <- "vehiclesbyroute"
  base_url <- paste("http://realtime.mbta.com/developer/api/v2/",query,"?api_key=",api_key,"&format=json",sep="")
  full_url <- paste(base_url,"&route=",route_id,sep="")
  rawdata <- readLines(full_url, warn = F)
  dl <- jsonlite::fromJSON(txt=rawdata,simplifyDataFrame = T)
  allout <- NULL
  for(i in 1:length(dl$direction$direction_id)){
    thisout <- data.frame(
      route_id=dl$route_id,
      route_name=dl$route_name,
      route_type=dl$route_type,
      mode_name=dl$mode_name,
      direction_id=dl$direction$direction_id[i],
      direction_name=dl$direction$direction_name[i],
      trip_id=dl$direction$trip[[i]]$trip_id,
      trip_name=dl$direction$trip[[i]]$trip_name,
      trip_headsign=dl$direction$trip[[i]]$trip_headsign,
      vehicle_id=dl$direction$trip[[i]]$vehicle$vehicle_id,
      vehicle_lat=dl$direction$trip[[i]]$vehicle$vehicle_lat,
      vehicle_lon=dl$direction$trip[[i]]$vehicle$vehicle_lon,
      vehicle_bearing=ifelse(length(dl$direction$trip[[i]]$vehicle$vehicle_bearing)>0,dl$direction$trip[[i]]$vehicle$vehicle_bearing,NA), # not available for bus
      vehicle_timestamp=dl$direction$trip[[i]]$vehicle$vehicle_timestamp)
    allout <- rbind(allout,thisout)
  }
  return(allout)
}
