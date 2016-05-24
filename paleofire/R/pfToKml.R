#' Export selected site to Google Earth kml format
#' 
#' Export sites selected using pfSiteSel function to Google Earth kml format.
#' 
#' 
#' @param x An object of the class "pfSiteSel"
#' @param file File location and name with kml extension e.g.
#' file="/Users/Olivier/Desktop/truc.kml"
#' @return No value returned.
#' @author O. Blarquez
#' @examples
#' 
#' \dontrun{
#' x=pfSiteSel(id_site==222)
#' pfToKml(x,file="site222.kml")
#' }
#' 
pfToKml=function(x,file="NULL")
{
  if(file=="NULL") stop("Output not specified")

## For testing
# x=pfSiteSel(lat> 40, long>-85, long<(-70), lat<60)
# plot(x) ...
# file="/Users/Olivier/Desktop/truc.kml"

df=data.frame(name=row.names(summary(x)),summary(x))
sp::coordinates(df) = ~long+lat

sp::proj4string(df)<-sp::CRS("+init=epsg:4326")

options(warn=-1)
rgdal::writeOGR(df, dsn=file, layer= "df", driver="KML", dataset_options=c("NameField=name"))

}
