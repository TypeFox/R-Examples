Troutes <- function(api_key){
  # returns list of MBTA route ids and route names
  query <- "routes"
  base_url <- paste("http://realtime.mbta.com/developer/api/v2/",query,"?api_key=",api_key,"&format=json",sep="")
  full_url <- paste(base_url)
  rawdata <- readLines(full_url, warn = F)
  dl <- jsonlite::fromJSON(txt=rawdata,simplifyDataFrame = T,flatten=F)
  allout <- NULL
  for(i in 1:length(dl$mode$route)){
    route_types <- dl$mode$route_type[i]
    mode_names <- dl$mode$mode_name[i]
    routes <- dl$mode$route[[i]][,c("route_id","route_name")]
    thisout <- data.frame(route_type=route_types,mode_name=mode_names,routes)
    allout <- rbind(allout,thisout)
  }
  return(allout)
}

