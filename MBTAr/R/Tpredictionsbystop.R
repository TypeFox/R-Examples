Tpredictionsbystop <- function(stop_id,include_access_alerts=FALSE,include_service_alerts=TRUE,api_key){
  # finds predicted arrivals and departures in the next hour for a particular stop
  query <- "predictionsbystop"
  base_url <- paste("http://realtime.mbta.com/developer/api/v2/",query,"?api_key=",api_key,"&format=json",sep="")
  full_url <- paste(base_url,"&stop=",stop_id,"&include_access_alerts=",include_access_alerts,"&include_service_alerts=",include_service_alerts,sep="")
  rawdata <- readLines(full_url, warn = F)
  dl <- jsonlite::fromJSON(txt=rawdata,simplifyDataFrame = T)
  allout <- NULL
  for(j in 1:length(dl$mode$route[[1]]$route_id)){
    thisout <- NULL
    for(i in 1:length(dl$mode$route[[1]]$direction[[j]]$trip)){
      if(length(dl$mode$route[[1]]$direction[[j]]$trip[[i]]$vehicle)>0){
        thisroute <- data.frame(stop_id=dl$stop_id,
                                stop_name=dl$stop_name,
                                route_type=dl$mode$route_type,
                                mode_name=dl$mode$mode_name,
                                route_id=dl$mode$route[[1]]$route_id[j],
                                route_name=dl$mode$route[[1]]$route_name[j],
                                direction_id=dl$mode$route[[1]]$direction[[j]]$direction_id[i],
                                direction_name=dl$mode$route[[1]]$direction[[j]]$direction_name[i],
                                trip_id=dl$mode$route[[1]]$direction[[j]]$trip[[i]]$trip_id,
                                trip_name=dl$mode$route[[1]]$direction[[j]]$trip[[i]]$trip_name,
                                trip_headsign=dl$mode$route[[1]]$direction[[j]]$trip[[i]]$trip_headsign,
                                sch_arr_dt=dl$mode$route[[1]]$direction[[j]]$trip[[i]]$sch_arr_dt,
                                sch_dep_dt=dl$mode$route[[1]]$direction[[j]]$trip[[i]]$sch_dep_dt,
                                pre_dt=dl$mode$route[[1]]$direction[[j]]$trip[[i]]$pre_dt,
                                pre_away=dl$mode$route[[1]]$direction[[j]]$trip[[i]]$pre_away,
                                vehicle_id=dl$mode$route[[1]]$direction[[j]]$trip[[i]]$vehicle$vehicle_id,
                                vehicle_lat=dl$mode$route[[1]]$direction[[j]]$trip[[i]]$vehicle$vehicle_lat,
                                vehicle_lon=dl$mode$route[[1]]$direction[[j]]$trip[[i]]$vehicle$vehicle_lon,
                                vehicle_bearing=dl$mode$route[[1]]$direction[[j]]$trip[[i]]$vehicle$vehicle_bearing,
                                vehicle_timestamp=dl$mode$route[[1]]$direction[[j]]$trip[[i]]$vehicle$vehicle_timestamp)
      }
      if(length(dl$mode$route[[1]]$direction[[j]]$trip[[i]]$vehicle)==0){
        thisroute <- data.frame(stop_id=dl$stop_id,
                                stop_name=dl$stop_name,
                                route_type=dl$mode$route_type,
                                mode_name=dl$mode$mode_name,
                                route_id=dl$mode$route[[1]]$route_id[j],
                                route_name=dl$mode$route[[1]]$route_name[j],
                                direction_id=dl$mode$route[[1]]$direction[[j]]$direction_id[i],
                                direction_name=dl$mode$route[[1]]$direction[[j]]$direction_name[i],
                                trip_id=dl$mode$route[[1]]$direction[[j]]$trip[[i]]$trip_id,
                                trip_name=dl$mode$route[[1]]$direction[[j]]$trip[[i]]$trip_name,
                                trip_headsign=dl$mode$route[[1]]$direction[[j]]$trip[[i]]$trip_headsign,
                                sch_arr_dt=dl$mode$route[[1]]$direction[[j]]$trip[[i]]$sch_arr_dt,
                                sch_dep_dt=dl$mode$route[[1]]$direction[[j]]$trip[[i]]$sch_dep_dt,
                                pre_dt=dl$mode$route[[1]]$direction[[j]]$trip[[i]]$pre_dt,
                                pre_away=dl$mode$route[[1]]$direction[[j]]$trip[[i]]$pre_away,
                                vehicle_id=NA,
                                vehicle_lat=NA,
                                vehicle_lon=NA,
                                vehicle_bearing=NA,
                                vehicle_timestamp=NA)
      }
      thisout <- rbind(thisout,thisroute)
    }
    allout <- rbind(allout,thisout)
  }
  return(allout)
}
