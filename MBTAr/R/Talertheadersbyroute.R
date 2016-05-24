Talertheadersbyroute <- function(route_id,include_access_alerts=FALSE,include_service_alerts=TRUE,api_key){
  query <- "alertheadersbyroute"
  base_url <- paste("http://realtime.mbta.com/developer/api/v2/",query,"?api_key=",api_key,"&format=json",sep="")
  full_url <- paste(base_url,"&route=",route_id,"&include_access_alerts=",include_access_alerts,"&include_service_alerts=",include_service_alerts,sep="")
  rawdata <- readLines(full_url, warn = F)
  dl <- jsonlite::fromJSON(txt=rawdata,simplifyDataFrame = T,flatten=F)
  allout <- data.frame(route_id = dl$route_id,
                       route_name = dl$route_name,
                       alert_id = dl$alert_headers$alert_id,
                       header_text = dl$alert_headers$header_text)
  return(allout)
}
