list_issue_types <- function(city = NULL, lat = NULL, long = NULL, limit = 100) {
  total <- 0
  page <- 1
  if(length(city)>0 & (length(lat)>0 | length(long)>0)){
    lat <- NULL
    long <- NULL
    warning("Cannot specify both city and lat/long locations. Using city...")
  }
  if((length(lat)>0 & length(long)<1) | length(lat)<1 & length(long)>0){
    stop("Specify valid lat/long pair or city")
  }
  url <- paste("https://seeclickfix.com/api/v2/issues/new?",ifelse(length(city)>0,paste("address=",city,sep=""),""),ifelse(length(lat)>0,paste("lat=", lat,"&lng=",long,sep=""),""),"&per_page=",limit,"&page=",page, sep = "")
  url <- gsub(" ","%20",x=url)
  rawdata <- RCurl::getURL(url)
  scf <- jsonlite::fromJSON(txt=rawdata,simplifyDataFrame = T,flatten=F)
  allout <- data.frame(title = scf$request_types$title,
                       organization = scf$request_types$organization,
                       url = scf$request_types$url,
                       potential_duplicate_issues_url = scf$request_types$potential_duplicate_issues_url
  )
  
  total <- nrow(allout)
  
  while(limit>total){
    page <- page+1
    if((limit-total)<100){limit <- (limit-total)}
    url <- paste("https://seeclickfix.com/api/v2/issues/new?address=", city,"&per_page=",limit,"&page=",page, sep = "")
    url <- gsub(" ","%20",x=url)
    rawdata <- RCurl::getURL(url)
    scf <- jsonlite::fromJSON(txt=rawdata,simplifyDataFrame = T,flatten=F)
    holder <- data.frame(title = scf$request_types$title,
                         organization = scf$request_types$organization,
                         url = scf$request_types$url,
                         potential_duplicate_issues_url = scf$request_types$potential_duplicate_issues_url
    )
    allout <- rbind(allout,holder)
    total <- nrow(allout) 
  }
  return(allout)
}
