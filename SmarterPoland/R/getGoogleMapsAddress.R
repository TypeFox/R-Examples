getGoogleMapsAddress <-  function(street = "Banacha 2", city = "Warszawa", country="Poland", positionOnly = TRUE, delay=1) {
  apiHttps  <- paste0("http://maps.googleapis.com/maps/api/geocode/json?address=", street, ",+",city,",+",country,"&sensor=true")
  res <- sapply(apiHttps, function(apiHttp) {
    Sys.sleep(delay)
    getGoogleMapsSignleAddress(apiHttp, positionOnly)
  })
  if (class(res) == "matrix") res <- t(res)
  res
}
  
getGoogleMapsSignleAddress <-
	function(apiHttp, positionOnly = TRUE) {
    level <- 0
    apiHttp <- gsub(apiHttp, pattern=" ", replacement="\\+")
    jsnip <-rjson::fromJSON( file=apiHttp, method = "C" )
    
    if (length(jsnip[[1]]) == 0) {
      apiHttp <- gsub(apiHttp, pattern="[0-9]", replacement="")
      jsnip <- rjson::fromJSON( file=apiHttp, method = "C" )
      level <- 1
      if (length(jsnip[[1]]) == 0) {
        apiHttp <- gsub(apiHttp, pattern="address=[^,]*,", replacement="address=")
        jsnip <- rjson::fromJSON( file=apiHttp, method = "C" )
        level <- 2
        if (length(jsnip[[1]]) == 0) {
          apiHttp <- gsub(apiHttp, pattern="address=[^,]*,", replacement="address=")
          jsnip <- rjson::fromJSON( file=apiHttp, method = "C" )
          level <- 3
        } 
      } 
    } 
  	if (length(jsnip) == 2 & jsnip$status == "OVER_QUERY_LIMIT") 
        return("OVER_QUERY_LIMIT")
    if (positionOnly)
       return(c(unlist(jsnip[[1]][[1]]$geometry$location),level))
	  jsnip$level = level
	  jsnip
}
