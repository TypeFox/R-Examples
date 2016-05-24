Tschedulebytrip <- function(trip_id,datetime=Sys.time(),api_key){
  # finds scheduled arrivals and departures for a given trip
  query <- "schedulebytrip"
  datetime <- as.integer(datetime)
  base_url <- paste("http://realtime.mbta.com/developer/api/v2/",query,"?api_key=",api_key,"&format=json",sep="")
  full_url <- paste(base_url,"&trip=",trip_id,"&datetime=",datetime,sep="")
  rawdata <- readLines(full_url, warn = F)
  dl <- jsonlite::fromJSON(txt=rawdata,simplifyDataFrame = T)
  allout <- data.frame(route_id=dl$route_id,
                       route_name=dl$route_name,
                       trip_id=dl$trip_id,
                       trip_name=dl$trip_name,
                       direction_id=dl$direction_id,
                       direction_name=dl$direction_name,
                       stop_sequence=dl$stop$stop_sequence,
                       stop_id=dl$stop$stop_id,
                       stop_name=dl$stop$stop_name,
                       sch_arr_dt=dl$stop$sch_arr_dt,
                       sch_dep_dt=dl$stop$sch_dep_dt)
  return(allout)
}
