list_places <- function(startingpoint,limit=100){
	total <- 0
	page <- 1
	url <- paste("https://seeclickfix.com/api/v2/places?address=", startingpoint, "&per_page=",limit,"&page=",page, sep = "")
	url <- gsub(" ","%20",x=url)
	rawdata <- RCurl::getURL(url)
  scf <- jsonlite::fromJSON(txt=rawdata,simplifyDataFrame = T,flatten=F)
  allout <- NULL
  for(i in 1:length(scf$places$center$coordinates)){
    thisout <- data.frame(id = scf$places$id[i],
                          name = scf$places$name[i],
                          url_name = scf$places$url_name[i],
                          county = scf$places$county[i],
                          state = scf$places$state[i],
                          place_type = scf$places$place_type[i],
                          url = scf$places$url[i],
                          html_url = scf$places$html_url[i],
                          html_report_url = scf$places$html_report_url[i],
                          type = scf$places$center$type[i],
                          lat = scf$places$center$coordinates[[i]][2],
                          lng = scf$places$center$coordinates[[i]][1]
    )
    allout <- rbind(allout,thisout)
  }
  total <- nrow(allout)

	while(limit>total){
		page <- page+1
		if((limit-total)<100){limit <- (limit-total)}
  url <- paste("https://seeclickfix.com/api/v2/places?address=", startingpoint, "&per_page=",(limit-total),"&page=",page, sep = "")
  url <- gsub(" ","%20",x=url)
  rawdata <- RCurl::getURL(url)
  scf <- jsonlite::fromJSON(txt=rawdata,simplifyDataFrame = T,flatten=F)
  holder <- NULL
  for(i in 1:length(scf$places$center$coordinates)){
    thisout <- data.frame(id = scf$places$id[i],
                          name = scf$places$name[i],
                          url_name = scf$places$url_name[i],
                          county = scf$places$county[i],
                          state = scf$places$state[i],
                          place_type = scf$places$place_type[i],
                          url = scf$places$url[i],
                          html_url = scf$places$html_url[i],
                          html_report_url = scf$places$html_report_url[i],
                          type = scf$places$center$type[i],
                          lat = scf$places$center$coordinates[[i]][2],
                          lng = scf$places$center$coordinates[[i]][1]
    )
    holder <- rbind(holder,thisout)
  }
  allout <- rbind(allout,holder)
  total <- nrow(allout)
  }
  return(allout)
}