#' Geodode an address by Google Maps
#'
#' @param gcStr A string
#' @return Geocoded address
#' @export
getGeoCodeCS <- function(gcStr)
{
  #OrgEnd<-gcStr
  gcStr <- iconv(gcStr,"latin1","UTF-8")
  gcStr <- gsub(' ','%20',gcStr) #Encode URL Parameters
  #Open Connection
  connectStr <- paste('http://maps.google.com/maps/api/geocode/json?sensor=false&language=pt-BR&address=',gcStr, sep="")
  con <- url(connectStr)
  data.json <- fromJSON(paste(readLines(con), collapse=""))
  close(con)
  #Flatten the received JSON
  data.json <- unlist(data.json)
  lat <- as.double(data.json["results.geometry.location.lat"])
  lng <- as.double(data.json["results.geometry.location.lng"])
  Status <- data.json["status"]
  FormattedAddress <- iconv(as.character(data.json["results.formatted_address"]),"UTF-8","latin1")
  gcodes <- c(lat, lng, FormattedAddress, Status)
  names(gcodes) <- c("Lat", "Lng","Ender","Status")
  return(gcodes)
}
