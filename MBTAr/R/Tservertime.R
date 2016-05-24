Tservertime <- function(api_key){
  # returns MBTA API server time
  query <- "servertime"
  base_url <- paste("http://realtime.mbta.com/developer/api/v2/",query,"?api_key=",api_key,"&format=json",sep="")
  full_url <- paste(base_url)
  rawdata <- readLines(full_url, warn = F)
  dl <- jsonlite::fromJSON(txt=rawdata,simplifyDataFrame = T,flatten=F)
  allout <- data.frame(server_dt=dl$server_dt)
  return(allout)
}

