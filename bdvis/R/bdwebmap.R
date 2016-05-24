#' Interactive web page based map of records
#' 
#' Generates a interactive web map and opens it in a new browser window 
#' 
#' @import leafletR
#' @importFrom utils browseURL
#' @param indf input data frame containing biodiversity data set
#' @examples \dontrun{
#'  bdwebmap(inat)
#' }
#' @family Spatial visualizations
#' @export
bdwebmap <- function(indf){
bs=indf[,c(which(names(indf)=="Scientific_name"),
           which(names(indf)=="Latitude"),
           which(names(indf)=="Longitude"))]
bs=bs[which(bs$Latitude>0),]
bs$Latitude=as.numeric(bs$Latitude)
bs$Longitude=as.numeric(bs$Longitude)
bsj=toGeoJSON(bs, lat.lon=c(2,3))

cols <- c("#8D5C00", "#2F5CD7","#E91974", "#3CB619","#7EAFCC",
          "#4F2755","#F5450E","#264C44","#3EA262","#FA43C9","#6E8604","#631D0E","#EE8099","#E5B25A",
          "#0C3D8A","#9E4CD3","#195C7B","#9F8450","#7A0666","#BBA3C5","#F064B4","#108223","#553502",
          "#17ADE7","#83C445","#C52817","#626302","#9F9215","#6CCD78","#BF3704")
pal <- cols[1:5]
sty <- styleCat(prop = "Scientific_name", val = unique(bs$Scientific_name), style.val = pal, leg = "Scientific Name")
unlink(".\\bs\\bs.html")
unlink(".\\bs\\bs.geojson")

map <- leaflet(bsj,base.map="tls",style=sty,popup = "Scientific_name")
browseURL(map)
}
