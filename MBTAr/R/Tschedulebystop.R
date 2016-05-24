Tschedulebystop <- function(stop_id,route_id,direction=NULL,datetime=Sys.time(),max_time=60,max_trips=5,api_key){
  # finds scheduled trips within a time window at a specified stop on a specified route
  query <- "schedulebystop"
  datetime <- as.integer(datetime)
  base_url <- paste("http://realtime.mbta.com/developer/api/v2/",query,"?api_key=",api_key,"&format=json",sep="")
  full_url <- paste(base_url,"&stop=",stop_id,ifelse(length(route_id)>0,paste("&route=",route_id,sep=""),""),ifelse(length(direction)>0,paste("&direction=",direction,sep=""),""),"&datetime=",datetime,"&max_time=",max_time,"&max_trips=",max_trips,sep="")
  rawdata <- readLines(full_url, warn = F)
  dl <- jsonlite::fromJSON(txt=rawdata,simplifyDataFrame = T)
  allout <- NULL
  for(j in 1:length(dl$mode$route[[1]]$route_id)){
    thisroute <- NULL
    for(k in 1:length(dl$mode$route[[1]]$direction[[j]]$direction_id)){
    stop_id=dl$stop_id
    stop_name=dl$stop_name
    if(length(dl$mode)<1){
      stop("Invalid stop_id/route_id/direction combination.\nTry direction=NULL or check Troutesbystop() to make sure that you have specified a combination that exists.")
    }
    route_type=dl$mode$route_type
    mode_name=dl$mode$mode_name
    route_id=dl$mode$route[[1]]$route_id[j]
    route_name=dl$mode$route[[1]]$route_name[j]
    direction_id=dl$mode$route[[1]]$direction[[j]]$direction_id[k]
    direction_name=dl$mode$route[[1]]$direction[[j]]$direction_name[k]
    trip_id=dl$mode$route[[1]]$direction[[j]]$trip[[k]]$trip_id
    trip_name=dl$mode$route[[1]]$direction[[j]]$trip[[k]]$trip_name
    sch_arr_dt=dl$mode$route[[1]]$direction[[j]]$trip[[k]]$sch_arr_dt
    sch_dep_dt=dl$mode$route[[1]]$direction[[j]]$trip[[k]]$sch_dep_dt
    thisdirection <- data.frame(stop_id,stop_name,route_type,mode_name,route_id,route_name,direction_id,direction_name,trip_id,trip_name,sch_arr_dt,sch_dep_dt)
    thisroute <- rbind(thisroute,thisdirection)
  }
    allout <- rbind(allout,thisroute)
  }
  return(allout)
}
