#' aoristic graph by shapefile boundary
#' @param spdf spatial point data frame produced from aoristic.spdf
#' @param area.shp spatial polygon data frame used as a boundary in WGS84
#' @param output a character representing the name of an output folder (the folder will be generated in the current working directory; default name = output) 
#' @return kml file
#' @references Ratcliffe, J. H. (2002). Aoristic Signatures and the Spatio-Temporal Analysis of High Volume Crime Patterns. Journal of Quantitative Criminology, 18(1), 23-43. 
#' @import lubridate classInt reshape2 GISTools ggplot2 spatstat plotKML RColorBrewer 
#' @importFrom sp proj4string
#' @export
#' @examples
#' \donttest{
#' data(aoristic)
#' data.spdf <- aoristic.spdf(data=arlington, 
#'    DateTimeFrom="DateTimeFrom", DateTimeTo="DateTimeTo", 
#'    lon="lon", lat="lat")
#' aoristic.shp(spdf=data.spdf, area.shp=CouncilDistrict)
#' }

aoristic.shp <- function(spdf, area.shp, output="output"){
  
  #defining variables (to avoid "Note" in the package creation)
  sortID=NULL
  time0=NULL
  time23=NULL
  freq=NULL
  
  # check projections
  m <- match.call()
  if (is.na(area.shp@proj4string@projargs) | area.shp@proj4string@projargs==""){
    stop(paste("Projection information of", m$area.shp, "is missing!"))
  }
  
  # ver 0.3
  CRS <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  # if(!area.shp@proj4string@projargs==CRS(CRS)@projargs){stop("the coordinate reference system is not in WGS84")}
  
  # ver 0.5
  # CRS <- "+init=epsg:4326"
  # if(!area.shp@proj4string@projargs==CRS(CRS)@projargs){
  #  area.shp <- spTransform(area.shp, CRS(CRS))
  #}
  
  if (!check_projection(area.shp)){
    area.shp <- reproject(area.shp)
  }
  # projection needs to be identical in words for spatial aggregation to work
  if (!area.shp@proj4string@projargs==spdf@proj4string@projargs){
    area.shp <- suppressMessages(reproject(area.shp, spdf@proj4string@projargs, show.output.on.console = FALSE))
  }
  
  # create output location
  folder.location <- getwd()
  tryCatch({
    dir.create(file.path(folder.location, output), showWarnings = FALSE)
    dir.create(file.path(folder.location, output, "GISboundary"), showWarnings = TRUE)
  }, warning=function(w){
    stop(paste("The output folder already exists in: ", getwd(), "/", output, "/GISboundary", sep=""))
  })
  
  # set output location
  setwd(file.path(folder.location, output, "GISboundary"))
  
  area.shp@data$sortID <- seq(1, nrow(area.shp@data), 1)
  area.shp <- area.shp[order(area.shp@data$sortID),]
  
  
  # aggregate aoristic count per GIS boundary through for-loop----------
  for (i in 4:27){
    agg <- aggregate(spdf[i], area.shp, FUN=sum)
    agg@data[is.na(agg@data)] <- 0
    area.shp <- spCbind(area.shp, agg@data)
  }
  
  # graph by areas ------------
  graph2 <- subset(area.shp@data, select=c(sortID, time0:time23))
  
  graph2 <- melt(graph2, id.vars="sortID")
  graph2$variable <- as.numeric(gsub("time", "", graph2$variable))
  # tail(graph2)
  
  # factor 
  graph2$variable <- factor(graph2$variable, 
                            levels=c("6", "7", "8", "9", "10", "11", "12", 
                                     "13", "14", "15", "16", "17", "18", "19", "20",
                                     "21", "22", "23", "0", "1", "2", "3", "4", "5"))
  
  names(graph2)[2:3] <- c("hour", "freq")
  # sort by area and hour
  graph2 <- graph2[order(graph2[1], graph2$hour),]
  
  # set output directory
  setwd(file.path(folder.location, output, "GISboundary"))
  
  # aoristic graphs by GIS area through for-loop
  for (i in 1:nrow(area.shp@data)){
    
    graph.temp<-graph2[graph2$sortID==i,]
    p <- ggplot(graph.temp, aes(x=hour, y=freq)) + 
      geom_bar(stat="identity") +
      ylim(0, max(graph2$freq))
    
    ggsave(filename=paste("sortID_", i, ".png", sep=""), plot=p, width = 6, height = 4)
    area.shp@data$img[i] <-  paste("sortID_", i, ".png", sep="")
  }
  
  ####
  ## create KML (shapefile)--------
  
  if (nrow(area.shp@data)>9){
    nclr <- 9 # the number of classification categories
  } else {
    nclr <- round(nrow(area.shp@data)/2)
  }
  
  plotclr <- brewer.pal(nclr,"YlOrRd")
  
  area.shp@data$Total <- rowSums(area.shp@data[,c("time0", "time1", "time2", "time3", "time4", 
                                                  "time5", "time6", "time7", "time8", "time9", 
                                                  "time10", "time11", "time12", "time13", "time14",
                                                  "time15", "time16", "time17", "time18", "time19",
                                                  "time20", "time21", "time22", "time23")])
  plotvar <- area.shp@data$Total
  
  class <- classIntervals(plotvar, nclr, style="jenks") 
  colcode <- findColours(class, plotclr, digits=4)
  area.shp@data$col <- add.alpha(colcode, 0.75)
  
  # relative path                          
  out <- sapply(slot(area.shp, "polygons"), function(x) { kmlPolygon(x,
                                                                     # name=as(area.shp, "data.frame")[slot(x, "ID"), "DISTRICTDE"], 
                                                                     name=paste("Crime Count: ", round(as(area.shp, "data.frame")[slot(x, "ID"), "Total"]), sep=""), 
                                                                     col =as(area.shp, "data.frame")[slot(x, "ID"), "col"], lwd=1, border='black', 
                                                                     description=paste("<img src=", 
                                                                                       as(area.shp, "data.frame")[slot(x, "ID"), "img"], " width=\"450\">", sep=""))})
  
  kml.folder <- file.path(folder.location, output, "GISboundary")
  tf <- file.path(kml.folder, "Aoristic_GIS_boundary.kml")
  
  kmlFile <- file(tf, "w")
  cat(kmlPolygon(kmlname="Aoristic_GIS_boundary")$header,
      file=kmlFile, sep="\n")
  cat(unlist(out["style",]), file=kmlFile, sep="\n")
  cat(unlist(out["content",]), file=kmlFile, sep="\n")
  cat(kmlPolygon()$footer, file=kmlFile, sep="\n")
  close(kmlFile)
  
 setwd(folder.location) 
 
 print(paste("KML output file is in ", getwd(), "/", output, "/GISboundary", sep=""))
 
 browseURL(file.path(folder.location, output))
 
 browseURL(file.path(folder.location, output, "GISboundary", "Aoristic_GIS_boundary.kml"))
 
}