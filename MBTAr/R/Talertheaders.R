Talertheaders <- function(include_access_alerts=FALSE,include_service_alerts=TRUE,api_key){
  query <- "alertheaders"
  base_url <- paste("http://realtime.mbta.com/developer/api/v2/",query,"?api_key=",api_key,"&format=json",sep="")
  full_url <- paste(base_url,"&include_access_alerts=",include_access_alerts,"&include_service_alerts=",include_service_alerts,sep="")
  rawdata <- readLines(full_url, warn = F)
  dl <- jsonlite::fromJSON(txt=rawdata,simplifyDataFrame = T,flatten=F)
  allout <- data.frame(alert_id=dl$alert_headers$alert_id,header_text=dl$alert_headers$header_text)
  return(allout)
}
