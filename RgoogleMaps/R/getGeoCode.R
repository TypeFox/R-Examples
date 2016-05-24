getGeoCode <- structure(function#geocoding utility
### Geocode your data using, R, JSON and Google Maps' Geocoding APIs 
### see http://allthingsr.blogspot.de/2012/01/geocode-your-data-using-r-json-and.html
### and 
(
  gcStr, ##<< adddress to geocode
  verbose=0 ##<< level of verbosity
){
  #library("RJSONIO") #Load Library
  gcStr <- enc2utf8(gsub(' ','%20',gcStr)) #Encode URL Parameters
  #Open Connection
  connectStr <- paste('http://maps.google.com/maps/api/geocode/json?sensor=false&address=',gcStr, sep="") 
  if (verbose) cat("fetching ", connectStr, "\n")
  con <- url(connectStr)
  data.json <- fromJSON(paste(readLines(con), collapse=""))
  close(con)
  #Flatten the received JSON
  data.json <- unlist(data.json)
  lat <- data.json["results.geometry.location.lat"]
  lng <- data.json["results.geometry.location.lng"]
  gcodes <- as.numeric(c(lat, lng))
  names(gcodes) <- c("lat", "lon")
  return (gcodes)
### returns lat/lon for address
}, ex = function(){
  getGeoCode("Brooklyn")
  #You can run this on the entire column of a data frame or a data table:
  DF = cbind.data.frame(address=c("Berlin,Germany", "Princeton,NJ", 
            "cadillac+mountain+acadia+national+park"), lat = NA, lon = NA)
  DF <- with(DF, data.frame(address, t(sapply(DF$address, getGeoCode))))
  
})
