Tschedulebyroute <- function(route_id,direction=NULL,datetime=Sys.time(),max_time=60,max_trips=5,api_key){
  # finds scheduled trips within a time window on a specified route
  query <- "schedulebyroute"
  datetime <- as.integer(datetime)
  base_url <- paste("http://realtime.mbta.com/developer/api/v2/",query,"?api_key=",api_key,"&format=json",sep="")
  full_url <- paste(base_url,"&route=",route_id,ifelse(length(direction)>0,paste("&direction=",direction,sep=""),""),"&datetime=",datetime,"&max_time=",max_time,"&max_trips=",max_trips,sep="")
  rawdata <- readLines(full_url, warn = F)
  dl <- jsonlite::fromJSON(txt=rawdata,simplifyDataFrame = T)
  allout <- NULL
  for(k in 1:length(dl$direction$direction_id)){
    thisdirection <- NULL
      for(j in 1:length(dl$direction$trip[[k]]$trip_id)){
      route_id=dl$route_id
      route_name=dl$route_name
      direction_id=dl$direction$direction_id[k]
      direction_name=dl$direction$direction_name[k]
      trip_id=dl$direction$trip[[k]]$trip_id[j]
      trip_name=dl$direction$trip[[k]]$trip_name[j]
      stop_sequence=dl$direction$trip[[k]]$stop[[j]]$stop_sequence
      stop_id=dl$direction$trip[[k]]$stop[[j]]$stop_id
      stop_name=dl$direction$trip[[k]]$stop[[j]]$stop_name
      sch_arr_dt=stop_sequence=dl$direction$trip[[k]]$stop[[j]]$sch_arr_dt
      sch_dep_dt=stop_sequence=dl$direction$trip[[k]]$stop[[j]]$sch_dep_dt
      thistrip <- data.frame(route_id,route_name,direction_id,direction_name,trip_id,trip_name,stop_sequence,sch_arr_dt,sch_dep_dt)
      thisdirection <- rbind(thisdirection,thistrip)
    }
    allout <- rbind(allout,thisdirection)
  }
  return(allout)
}
