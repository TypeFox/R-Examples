#' Geocode addresses by Google Maps
#'
#' @param datafr A dataframe
#' @param addresscol A vector
#' @return The submited dataframe plus columns Lat, Lng, Address (Google Maps style) and Status
#'
#' @export
setLatLong <- function(datafr,addresscol)
{
  lista<-sapply(addresscol,function(val){getGeoCodeCS(val)},USE.NAMES=F)
  DF<-as.data.frame(matrix(t(lista),ncol(lista),4),stringsAsFactors=F)
  names(DF)<-c("Lat","Lng","End","Status")
  DF$Lat=as.double(DF$Lat)
  DF$Lng=as.double(DF$Lng)
  datafr<-data.frame(datafr,DF)
  return(datafr)
}
