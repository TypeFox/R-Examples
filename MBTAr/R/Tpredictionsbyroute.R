Tpredictionsbyroute <- function(route_id,include_access_alerts=FALSE,include_service_alerts=TRUE,api_key){
  query <- "predictionsbyroute"
  base_url <- paste("http://realtime.mbta.com/developer/api/v2/",query,"?api_key=",api_key,"&format=json",sep="")
  full_url <- paste(base_url,"&route=",route_id,"&include_access_alerts=",include_access_alerts,"&include_service_alerts=",include_service_alerts,sep="")
  rawdata <- readLines(full_url, warn = F)
  dl <- jsonlite::fromJSON(txt=rawdata,simplifyDataFrame = T)
  if(length(grep("Green",x=route_id))>0){
    warning("Requested predictions for Green Line - predictions and real time tracking not yet available for the Green Line.")
  }
  allout <- NULL
  for(i in 1:length(dl$direction$trip)){
    thisdirection <- NULL
    for(k in 1:length(dl$direction$trip[[i]]$trip_id)){
      thistrip <- data.frame(route_id=dl$route_id,
                             route_name=dl$route_name,
                             route_type=dl$route_type,
                             mode_name=dl$mode_name,
                             direction_id=dl$direction$direction_id[i],
                             direction_name=dl$direction$direction_name[i],
                             trip_id=dl$direction$trip[[i]]$trip_id[k],
                             trip_name=dl$direction$trip[[i]]$trip_name[k],
                             trip_headsign=dl$direction$trip[[i]]$trip_headsign[k],
                             vehicle_id=dl$direction$trip[[i]]$vehicle$vehicle_id[k],
                             vehicle_lat=dl$direction$trip[[i]]$vehicle$vehicle_lat[k],
                             vehicle_lon=dl$direction$trip[[i]]$vehicle$vehicle_lon[k],
                             vehicle_bearing=ifelse(length(dl$direction$trip[[i]]$vehicle$vehicle_bearing[k])>0,dl$direction$trip[[i]]$vehicle$vehicle_bearing[k],NA), # bus doesnt have this
                             vehicle_timestamp=dl$direction$trip[[i]]$vehicle$vehicle_timestamp[k],
                             # green line doesn't have tracking or schedule yet, so predictions are empty:
                             stop_sequence=ifelse(dl$route_type!=0,dl$direction$trip[[i]]$stop[[k]]$stop_sequence,NA),
                             stop_id=ifelse(dl$route_type!=0,dl$direction$trip[[i]]$stop[[k]]$stop_id,NA),
                             stop_name=ifelse(dl$route_type!=0,dl$direction$trip[[i]]$stop[[k]]$stop_name,NA),
                             sch_arr_dt=ifelse(dl$route_type!=0,dl$direction$trip[[i]]$stop[[k]]$sch_arr_dt,NA),
                             sch_dep_dt=ifelse(dl$route_type!=0,dl$direction$trip[[i]]$stop[[k]]$sch_dep_dt,NA),
                             pre_dt=ifelse(dl$route_type!=0,dl$direction$trip[[i]]$stop[[k]]$pre_dt,NA),
                             pre_away=ifelse(dl$route_type!=0,dl$direction$trip[[i]]$stop[[k]]$pre_away,NA))

      thisdirection <- rbind(thisdirection,thistrip) # all trips in one direction
    }
    allout <- rbind(allout,thisdirection)
  }
  return(allout)
}
