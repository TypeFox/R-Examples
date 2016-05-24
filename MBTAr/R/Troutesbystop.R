Troutesbystop <- function(stop_id=NULL,api_key){
  query <- "routesbystop"
  base_url <- paste("http://realtime.mbta.com/developer/api/v2/",query,"?api_key=",api_key,"&format=json",sep="")
  full_url <- paste(base_url,"&stop=",stop_id,sep="")
  rawdata <- readLines(full_url, warn = F)
  dl <- jsonlite::fromJSON(txt=rawdata,simplifyDataFrame = T,flatten=T)
  allmodes <- NULL
  for(i in 1:length(dl$mode$route_type)){
    allroutes <- data.frame(stop_id=dl$stop_id,
                            stop_name=dl$stop_name,
                            route_type=dl$mode$route_type[i],
                            mode_name=dl$mode$mode_name[i],
                            route_id=dl$mode$route[[i]]$route_id,
                            route_name=dl$mode$route[[i]]$route_name
    )
    allmodes <- rbind(allmodes,allroutes)
  }
  return(allroutes)
}

