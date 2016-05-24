CreateNeosComm <- function(curlopts = list(httpheader = c(`Content-Type` = "text/xml", 'User-Agent' = "R"), port = 3332), curlhandle = getCurlHandle()){
  url <- "http://www.neos-server.org"
  if(!("httpheader" %in% names(curlopts))){
    stop("\nNo 'httpheader' list element has been specified in 'curlopts'.\n")
  }
  if(!("port" %in% names(curlopts))){
    stop("\nNo 'port' list element has been specified in 'curlopts'.\n")
  }
  if(!(class(curlhandle) == "CURLHandle")){
    stop("\nObject for 'curlhandle' must be of class 'CURLHandle'.\n")
  }
  result <- new("NeosComm", url = url, curlopts = curlopts, curlhandle = curlhandle)
  return(result)
}
