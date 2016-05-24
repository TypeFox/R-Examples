Tstopsbylocation <- function(lat,lon,api_key){
  # finds MBTA stops by latitude and longitude coordinates. returns up to 15 within a 1-mile radius
  query <- "stopsbylocation"
  base_url <- paste("http://realtime.mbta.com/developer/api/v2/",query,"?api_key=",api_key,"&format=json",sep="")
  full_url <- paste(base_url,"&lat=",lat,"&lon=",lon,sep="")
  rawdata <- readLines(full_url, warn = F)
  dl <- jsonlite::fromJSON(txt=rawdata,simplifyDataFrame = T)
  allout <- dl$stop
  #stop_id, stop_name, parent_station, parent_station_name, stop_lat, stop_lon, distance
  return(allout)
}
