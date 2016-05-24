## Andrew Heiss
## Date: 31 January 2013
## Include project-specific functions in this file

.onLoad <- function(libname, pkgname){
  if (is.null(getOption("RDSTK_api_base"))) {
    default_base <- "http://www.datasciencetoolkit.org"
    options("RDSTK_api_base"=default_base)
  }
}

street2coordinates <- function(address, session=getCurlHandle()) {
  api <- paste(getOption("RDSTK_api_base"), "/street2coordinates/", sep="")
  get.addy <- getURL(paste(api, URLencode(address), sep=""), curl=session)
  result <- ldply(fromJSON(get.addy), data.frame)
  names(result)[1] <- "full.address"
  return(result)
}

ip2coordinates <- function(ip, session=getCurlHandle()) {
  api <- paste(getOption("RDSTK_api_base"), "/ip2coordinates/", sep="")
  get.ips <- getURL(paste(api, URLencode(ip), sep=""), curl=session) 
  result <- ldply(fromJSON(get.ips), data.frame)
  names(result)[1] <- "ip.address"
  return(result)
}

coordinates2politics <- function(latitude, longitude, session=getCurlHandle()) {
  api <- paste(getOption("RDSTK_api_base"), "/coordinates2politics/", sep="")
  result <- getURL(paste(api, latitude, "%2c", longitude, sep=""), curl=session)
  return(result)
}

text2sentences <- function(text, session=getCurlHandle()) {
  api <- paste(getOption("RDSTK_api_base"), "/text2sentences/", sep="")
  r = dynCurlReader()
  curlPerform(postfields=text, url=api, post=1L, writefunction=r$update,
              curl=session)
  result <- fromJSON(r$value())
  return(result)
}

text2people <- function(text, session=getCurlHandle()) {
  api <- paste(getOption("RDSTK_api_base"), "/text2people/", sep="")
  r = dynCurlReader()
  curlPerform(postfields=text, url=api, post=1L, writefunction=r$update,
              curl=session)
  result <- ldply(fromJSON(r$value()), data.frame)
  return(result)
}

html2text <- function(html, session=getCurlHandle()) {
  api <- paste(getOption("RDSTK_api_base"), "/html2text/", sep="")
  r = dynCurlReader()
  curlPerform(postfields=html, url=api, post=1L, writefunction=r$update, 
              curl=session)
  result <- fromJSON(r$value())
  return(result)
}

text2times <- function(text, session=getCurlHandle()) {
  api <- paste(getOption("RDSTK_api_base"), "/text2times/", sep="")
  r = dynCurlReader()
  curlPerform(postfields=text, url=api, post=1L, writefunction=r$update,
              curl=session)
  result <- ldply(fromJSON(r$value()), data.frame)
  return(result)
}

text2sentiment <- function(text, session=getCurlHandle()) {
  api <- paste(getOption("RDSTK_api_base"), "/text2sentiment/", sep="")
  r = dynCurlReader()
  curlPerform(postfields=text, url=api, post=1L, writefunction=r$update,
              curl=session)
  result <- fromJSON(r$value())
  return(result)
}

coordinates2statistics <- function(latitude, longitude, statistic, session=getCurlHandle()) {
  api <- paste(getOption("RDSTK_api_base"), "/coordinates2statistics/", sep="")
  r <- getURL(paste(api, latitude, "%2c", longitude, "?statistics=", statistic, sep=""), curl=session)
  result <- ldply(fromJSON(r), data.frame)
  return(result)
}

