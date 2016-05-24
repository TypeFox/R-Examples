mixpanelUpdateProfile = function(
  account,  
  distinctID,  
  data, 
  updateLastSeen=FALSE,
  updateLocation=FALSE, 
  retryCount=100 
) {
  data[["$token"]] = jsonlite::unbox(account$token)
  data[["$distinct_id"]] = jsonlite::unbox(distinctID)
  data[["$ignore_time"]] = jsonlite::unbox(!updateLastSeen)
  if (!updateLocation)
    data[["ip"]] = jsonlite::unbox(0)     
  
  jsonData = jsonlite::toJSON(list(data), digits=3)
  
  cat(substr(jsonData, 1, 100), "...\n", sep="")
  
  data64 = RCurl::base64(jsonData, encode=TRUE, mode="character")
  endpoint = 'http://api.mixpanel.com/'
  url = paste(endpoint, "engage", sep="")   #?data=", data64, sep="")
  
  for (i in 1:retryCount) {
    try({
      res = RCurl::getForm(url, data=data64, style='HTTPGET')
      if (res == "1") # If OK, exit loop.
        break
      
      Sys.sleep(2) ## Retry after sleep.
      cat("retry", i, "")
    })
  }
  cat(res, "\n")
}
