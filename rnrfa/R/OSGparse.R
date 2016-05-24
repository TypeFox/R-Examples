#' Converts OS Grid Reference to BNG/WGS coordinates.
#'
#' @author Claudia Vitolo
#'
#' @description This function converts an Ordnance Survey (OS) grid reference to easting/northing or latitude/longitude coordinates.
#'
#' @param gridRefs This is a string (or a character vector) that contains the OS grid Reference.
#' @param CoordSystem By default, this is "BNG" which stands for British National Grids. The other option is to set CoordSystem = "WGS84", which returns latitude/longitude coordinates (more info can be found here https://www.epsg-registry.org/).
#'
#' @return vector made of two elements: the easting and northing (by default) or latitude and longitude coordinates.
#'
#' @examples
#' # single entry
#' OSGparse(gridRefs="TQ722213")
#'
#' # multiple entries
#' OSGparse(gridRefs=c("SN831869","SN829838"))
#'

OSGparse <- function(gridRefs, CoordSystem = "BNG" ) {

  # require(sp)

  # gridRefs <- as.list(gridRefs)
  xlon <- c()
  ylat <- c()

  for (gridRef in gridRefs){

    # For testinG: library(sp); library(rgdal); gridRef <- "SN949824"

    # Starting point is the south-west corner

    # First letter identifies the 500x500 km grid
    firstLetter <- substr(gridRef,1,1)
    # Englan + Wales + Scotland
    if (firstLetter=="S") {xOffset1 <- 0;   yOffset1 <- 0}
    if (firstLetter=="T") {xOffset1 <- 500; yOffset1 <- 0}
    if (firstLetter=="N") {xOffset1 <- 0;   yOffset1 <- 500}
    if (firstLetter=="H") {xOffset1 <- 0;   yOffset1 <- 1000}
    if (firstLetter=="O") {xOffset1 <- 500; yOffset1 <- 500}

    # Norther Ireland?
    if (firstLetter=="I") {
      yOffset1 <- 0
      xOffset1 <- 0
      EPSG <- 29902
      defaultCRS <- CRS("+init=epsg:29902")
    }else{
      EPSG <- 27700
      defaultCRS <- CRS("+init=epsg:27700")
    }

    # Second letter identifies the 100x100 km grid
    secondLetter <- substr(gridRef,2,2)
    if (secondLetter=="A") {yOffset2 <- 400; xOffset2 <- 0}
    if (secondLetter=="B") {yOffset2 <- 400; xOffset2 <- 100}
    if (secondLetter=="C") {yOffset2 <- 400; xOffset2 <- 200}
    if (secondLetter=="D") {yOffset2 <- 400; xOffset2 <- 300}
    if (secondLetter=="E") {yOffset2 <- 400; xOffset2 <- 400}
    if (secondLetter=="F") {yOffset2 <- 300; xOffset2 <- 0}
    if (secondLetter=="G") {yOffset2 <- 300; xOffset2 <- 100}
    if (secondLetter=="H") {yOffset2 <- 300; xOffset2 <- 200}
    if (secondLetter=="J") {yOffset2 <- 300; xOffset2 <- 300}
    if (secondLetter=="K") {yOffset2 <- 300; xOffset2 <- 400}
    if (secondLetter=="L") {yOffset2 <- 200; xOffset2 <- 0}
    if (secondLetter=="M") {yOffset2 <- 200; xOffset2 <- 100}
    if (secondLetter=="N") {yOffset2 <- 200; xOffset2 <- 200}
    if (secondLetter=="O") {yOffset2 <- 200; xOffset2 <- 300}
    if (secondLetter=="P") {yOffset2 <- 200; xOffset2 <- 400}
    if (secondLetter=="Q") {yOffset2 <- 100; xOffset2 <- 0}
    if (secondLetter=="R") {yOffset2 <- 100; xOffset2 <- 100}
    if (secondLetter=="S") {yOffset2 <- 100; xOffset2 <- 200}
    if (secondLetter=="T") {yOffset2 <- 100; xOffset2 <- 300}
    if (secondLetter=="U") {yOffset2 <- 100; xOffset2 <- 400}
    if (secondLetter=="V") {yOffset2 <- 0; xOffset2 <- 0}
    if (secondLetter=="W") {yOffset2 <- 0; xOffset2 <- 100}
    if (secondLetter=="X") {yOffset2 <- 0; xOffset2 <- 200}
    if (secondLetter=="Y") {yOffset2 <- 0; xOffset2 <- 300}
    if (secondLetter=="Z") {yOffset2 <- 0; xOffset2 <- 400}

    # The grid square gets you 100km resolution, the the first four numbers are x coordinate in 10m steps and the second four numbers are y coords in 10m steps. You have to lookup the grid square offset from the map of grid square. For NU, thats 4 along and 6 up. For your example, extract the coordinates thus:
    x=as.numeric(paste(substr(gridRef,3,5),0,sep=""))
    y=as.numeric(paste(substr(gridRef,6,8),0,sep=""))

    # create a spatial data frame with the coordinates. We multiply your four-figures to get metres, and add the X and Y grid offset (also in metres)
    xy = data.frame(x=(xOffset1+xOffset2)*1000 + x*10,
                    y=(yOffset1+yOffset2)*1000 + y*10)
    coordinates(xy) <- ~x+y
    proj4string(xy) <- defaultCRS

    # this is the epsg code for OSgrid references in metres. To convert to lat-long WGS84 coords:
    if (CoordSystem == "WGS84") {
      xy <- spTransform(xy,CRS("+init=epsg:4326"))
    }

    xlon <- c(xlon, xy@coords[[1]])
    ylat <- c(ylat, xy@coords[[2]])

  }

  if (CoordSystem == "BNG") newCoords <- list("easting"=xlon, "northing"=ylat)
  if (CoordSystem == "WGS84") newCoords <- list("lon"=xlon, "lat"=ylat)

  return(newCoords)

}
