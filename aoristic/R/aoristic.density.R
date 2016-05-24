#' aoristic graph by grid count
#' @param spdf spatial point data frame produced from aoristic.spdf 
#' @param h h parameter for the function kde2d (default=0.01); std. dev. 
#' @param n n parameter for the function kde2d (default=128); No. of cells in xy directions
#' @param probs percentile to be identified as hot spots for the function quantile (default=0.99)
#' @param output a character representing the name of an output folder (the folder will be generated in the current working directory; default name = output) 
#' @return kml file
#' @references Ratcliffe, J. H. (2002). Aoristic Signatures and the Spatio-Temporal Analysis of High Volume Crime Patterns. Journal of Quantitative Criminology, 18(1), 23-43. 
#' @import lubridate classInt reshape2 GISTools ggplot2 spatstat
#' @importFrom MASS kde2d
#' @importFrom sp proj4string Polygons Polygon SpatialPolygons CRS SpatialPolygonsDataFrame as.image.SpatialGridDataFrame
#' @importFrom maptools spCbind kmlPolygon GE_SpatialGrid kmlOverlay
#' @export
#' @examples
#' \donttest{
#' data(aoristic)
#' data.spdf <- aoristic.spdf(data=arlington, 
#'    DateTimeFrom="DateTimeFrom", DateTimeTo="DateTimeTo", 
#'    lon="lon", lat="lat")
#' aoristic.density(spdf=data.spdf)
#' }
aoristic.density <- function(spdf, h=0.01, n=128, probs=0.99, output="output"){
  
  #defining variables (to avoid "Note" in the package creation)
  sortID=NULL
  time0=NULL
  time23=NULL
  freq=NULL
  
  # creating output folder
  folder.location <- getwd()
  
  tryCatch({
    dir.create(file.path(folder.location, output), showWarnings = FALSE)
    dir.create(file.path(folder.location, output, "Density and Contour"), showWarnings = TRUE)
  }, warning=function(w){
    stop(paste("The output folder already exists in: ", getwd(), "/", output, "/Density and Contour", sep=""))
  })
  
  setwd(file.path(folder.location, output, "Density and Contour"))
  
  # create point data
  data.ppp <- as(spdf, "ppp")
  
  kde <- kde2d(x=data.ppp$x, y=data.ppp$y, h=h, n=n) 
  # image(kde)
  # quantile(kde$z, 0.99)
  # contour(kde, levels=c(quantile(kde$z, 0.99)), add=TRUE)
  # c <- contourLines(kde$x, kde$y, kde$z)
  c <- contourLines(kde, levels=c(quantile(kde$z, probs)))
  # convert a list to SpatialPolygons
  for (i in 1:length(c)){
    xy <- cbind(c[[i]]$x, c[[i]]$y)
    p <- Polygon(xy)
    ps <- Polygons(list(p), i)
    sps.temp <- SpatialPolygons(list(ps))
    
    if (i==1){
      c.sps <- sps.temp
    } else {
      c.sps <- rbind(c.sps, sps.temp) 
    }
  }
  sp::proj4string(c.sps) <- CRS(spdf@proj4string@projargs)
  
  id <- data.frame(id=seq(1:length(c.sps)))
  c.sps <- SpatialPolygonsDataFrame(c.sps, data=id)
  
  # create SPDF for kernel contour count ------    	
  area.shp <- c.sps
  area.shp@proj4string <- CRS(spdf@proj4string@projargs)
  area.shp <- as(area.shp, "SpatialPolygonsDataFrame")
  area.shp@data$dummy <- seq(1, length(area.shp), 1)
  names(area.shp@data) <- "sortID"
  
  # aggregate aoristic count through for-loop (Kernel Contour)----------
  
  for (i in 4:27){
    agg <- aggregate(spdf[i], area.shp, FUN=sum)
    agg@data[is.na(agg@data)] <- 0
    area.shp <- spCbind(area.shp, agg@data)
  }
  
  # graph by areas (Kernel Contour) ------------
  graph2 <- subset(area.shp@data, select=c(sortID, time0:time23))
  
  graph2 <- melt(graph2, id.vars="sortID")
  graph2$variable <- as.numeric(gsub("time", "", graph2$variable))
  
  # factor 
  graph2$variable <- factor(graph2$variable, 
                            levels=c("6", "7", "8", "9", "10", "11", "12", 
                                     "13", "14", "15", "16", "17", "18", "19", "20",
                                     "21", "22", "23", "0", "1", "2", "3", "4", "5"))
  
  names(graph2)[2:3] <- c("hour", "freq")
  # sort by area and hour
  graph2 <- graph2[order(graph2[1], graph2$hour),]
  
  for (i in 1:nrow(area.shp@data)){
    
    graph.temp<-graph2[graph2$sortID==i,]
    p <- ggplot(graph.temp, aes(x=hour, y=freq)) + 
      geom_bar(stat="identity") +
      ylim(0, max(graph2$freq))
    
    ggsave(filename=paste("sortID_", i, ".png", sep=""), plot=p, width = 6, height = 4)
    area.shp@data$img[i] <-  paste("sortID_", i, ".png", sep="")
  }
  
  ## create KML (Kernel Contour) --------
  
  area.shp@data$Total <- rowSums(area.shp@data[,c("time0", "time1", "time2", "time3", "time4", 
                                                  "time5", "time6", "time7", "time8", "time9", 
                                                  "time10", "time11", "time12", "time13", "time14",
                                                  "time15", "time16", "time17", "time18", "time19",
                                                  "time20", "time21", "time22", "time23")])
  
  out <- sapply(slot(area.shp, "polygons"), function(x) {kmlPolygon(x,
                                                                    name=paste("Crime Count: ", round(as(area.shp, "data.frame")[slot(x, "ID"), "Total"]), sep=""), 
                                                                    lwd=3, border='black', 
                                                                    description=paste("<img src=", 
                                                                                      as(area.shp, "data.frame")[slot(x, "ID"), "img"], " width=\"450\">", sep=""))})
  
  kml.folder <- file.path(folder.location, output, "Density and Contour")
  tf <- file.path(kml.folder, "Aoristic_Contour.kml")
  
  kmlFile <- file(tf, "w")
  
  cat(kmlPolygon(kmlname="Aoristic_Contour")$header,
      file=kmlFile, sep="\n")
  cat(unlist(out["style",]), file=kmlFile, sep="\n")
  cat(unlist(out["content",]), file=kmlFile, sep="\n")
  cat(kmlPolygon()$footer, file=kmlFile, sep="\n")
  close(kmlFile)
  
  
  # kernel density -> KML ------------------
  
  sp.pix <- kde.points(spdf, h=h, n=n)
  sp.grd <- as(sp.pix, "SpatialGridDataFrame")
  sp.grd@data$kde[sp.grd@data$kde < quantile(sp.grd@data$kde, 0.5)] <- NA
  
  #kernel2KML.R
  sp.grd.kml <- GE_SpatialGrid(sp.grd)
  
  tf <- file.path(kml.folder, "Density")
  
  png(filename=paste(tf, ".png", sep=""), width=sp.grd.kml$width,
      height=sp.grd.kml$height, bg="transparent")
  
  par(mar=c(0,0,0,0), xaxs="i", yaxs="i")
  
  Lab.palette <-
    colorRampPalette(c("green", "yellow", "red"), space = "Lab")
  
  image(as.image.SpatialGridDataFrame(sp.grd[1]), col=Lab.palette(10),
        xlim=sp.grd.kml$xlim, ylim=sp.grd.kml$ylim)
  
  kmlOverlay(sp.grd.kml, paste(tf, ".kml", sep=""), paste(tf, ".png", sep=""))
  dev.off()
  
  setwd(folder.location)
  
  print(paste("KML output file is in ", getwd(), "/", output, "/Density and Contour", sep=""))
  
  browseURL(file.path(folder.location, output))      
        
  browseURL(file.path(folder.location, output, "Density and Contour", "Density.kml"))
  Sys.sleep(10)
  browseURL(file.path(folder.location, output, "Density and Contour", "Aoristic_Contour.kml"))
        
  
}