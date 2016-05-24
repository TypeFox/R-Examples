#' Raw Data Class
#' 
#' A class to handle geoelectrics raw data.
#'
#' @slot address address of the raw data ascii file.
#' @slot seaLevel data frame that contains raw data positions and resitance values.
#' @slot height data frame that contains topography information (distance and height).
#' @export
#' @examples 
#' # rawData = new("RawData",
#' #                address = "../example/rawdata/p1_DipolDipol_SW-NE.dat")
#' 
#' data(sinkhole)
#' sinkhole@profiles[[2]]@rawData
#' sinkhole@profiles[[2]]@rawData@address
#' sinkhole@profiles[[2]]@rawData@height
#' sinkhole@profiles[[2]]@rawData@seaLevel
#' @seealso \code{\link{Profile-class}}, \code{\link{ProfileSet-class}}
setClass("RawData",
         representation = representation(
           address = "character",
           seaLevel = "data.frame",
           height = "data.frame"))
setMethod("initialize", "RawData",
          function(.Object, address) {
            if(nchar(address) == 0) {
              print("Raw Data address is missing.")
            } else {
              .Object@address = address
              con  <- file(address, open = "r")
              
              skipLines1 <- 9
              skipLines2 <- 0
              numberOfRows1 <- 0
              numberOfRows2 <- 0
              
              for(i in 1:10) {
                oneLine <- readLines(con, n=1)
              }
              
              while(grepl(".", oneLine, fixed=T)) {
                oneLine <- readLines(con, n=1)
                numberOfRows1 <- numberOfRows1 + 1
              }
              
              profile <- read.table(file=address, skip=skipLines1, 
                                    header=F, nrows=numberOfRows1)
              
              .Object@seaLevel <- data.frame(
                "dist"=profile[1],
                "depth"=profile[2],
                "val"=profile[4])
              colnames(.Object@seaLevel) <- c("dist", "depth", "val")
              
              skipLines2 <- skipLines1 + numberOfRows1
              
              try({
                while(!grepl(".", oneLine, fixed=T)) {
                  oneLine <- readLines(con, n=1)
                  skipLines2 <- skipLines2 + 1
                }
                
                while(grepl(".", oneLine, fixed=T)) {
                  oneLine <- readLines(con, n=1)
                  numberOfRows2 <- numberOfRows2 + 1
                }              
                
                height <- read.table(file=address, skip=skipLines2, 
                                     header=F, nrows=numberOfRows2)
                
                .Object@height <- data.frame(
                  "dist"=height[1], 
                  "height"=height[2])
                colnames(.Object@height) <- c("dist", "height")
              })
              
              close(con)              
            }
            return(.Object)
          })

#' XYZ Data Class
#' 
#' A class to handle xyz data. 
#' The software Res2DInv produces .xyz-files that contain the
#' inverted resistance values. The xyz class parses .xyz files.
#'
#' @slot address address of the xyz ascii file
#' @slot seaLevel data frame that contains positions and values withouth height adjustment
#' @slot heightAdaption data frame that contains positions and values after height adjustment 
#' @slot height data frame that contains topography information (distances and heights). 
#' It is reconstructed from .xyz-file.
#' @slot minData minimum value
#' @slot maxData maximum value
#' @export
#' @seealso \code{\link{Profile-class}}, \code{\link{ProfileSet-class}}, 
#' \code{\link{plotXyz}}, \code{\link{plotXyzHeight}}, \code{\link{plot3dXyz}}
#' @examples 
#' # xyzData = new("XyzData", 
#' #                address = "../example/xyzFiles/p1_DipolDipol_SW-NE.xyz"),
#' 
#' data(sinkhole) 
#' sinkhole@profiles[[1]]@xyzData
#' sinkhole@profiles[[1]]@xyzData@seaLevel
#' sinkhole@profiles[[1]]@xyzData@heightAdaption
#' sinkhole@profiles[[1]]@xyzData@height
#' sinkhole@profiles[[1]]@xyzData@minData
#' sinkhole@profiles[[1]]@xyzData@maxData
setClass("XyzData",
         representation = representation(
           address = "character",
           seaLevel = "data.frame",
           heightAdaption = "data.frame",
           minData = "numeric",
           maxData = "numeric",
           height = "data.frame"))
setMethod("initialize", "XyzData",
          function(.Object, address) {
            if(nchar(address) == 0) {
              print("XYZ data address is missing.")
            } else {
              .Object@address = address
              con  <- file(address, open = "r")
              
              skipLines1 <- 0
              numberOfRows <- 0
              numberOfRows2 <- 0
              
              oneLine <- readLines(con, n=1)
              while(grepl("/", oneLine)) {
                oneLine <- readLines(con, n=1)
                skipLines1 <- skipLines1 + 1
              }
              
              while(!grepl("/", oneLine)) {
                oneLine <- readLines(con, n=1)
                numberOfRows <- numberOfRows + 1    
              }
              
              skipLines2 <- skipLines1 + numberOfRows
              
              while(grepl("/", oneLine)) {
                oneLine <- readLines(con, n=1)
                skipLines2 <- skipLines2 + 1
              }
              
              while(!grepl("/", oneLine)) {
                oneLine <- readLines(con, n=1)
                numberOfRows2 <- numberOfRows2 + 1    
              }
              
              profile_without_topo <- read.table(file=address, skip=skipLines1, 
                                                 header=F, nrows=numberOfRows)
              .Object@seaLevel <- data.frame(
                dist=profile_without_topo[1],
                depth=profile_without_topo[2],
                val=profile_without_topo[3])
              colnames(.Object@seaLevel) <- c("dist", "depth", "val")
              
              profile <- read.table(file=address, skip=skipLines2, 
                                    header=F, nrows=numberOfRows2)
              
              .Object@heightAdaption <- data.frame(
                dist=profile[1],
                depth=profile[2],
                val=profile[3])
              colnames(.Object@heightAdaption) <- c("dist", "depth", "val")
              
              .Object@minData <- min(profile[3])
              .Object@maxData <- max(profile[3])
              
              height <- data.frame(
                dist=1,
                height=1)
            
              j <- 1
              for (i in 1:max(profile[1])) {
                indices <- which(round(profile[1])==i)
                if (length(indices)>0) {
                  index <- min(indices)
                  height[j,] <- c(profile[index,1], profile[index,2])
                  j <- j + 1
                }
              }
              
              .Object@height <- height
              
              close(con)
            }
            return(.Object)               
          })

#' GPS Coordinates Class
#'
#' A class to handle gps coordinates.
#'
#' @slot address address of the gps ascii file
#' @slot exact data frame that contains measured gps coordinates
#' @slot lm linear model of the measured gps coordinates
#' @slot relative relative coordinates
#' @slot lmRelative linear model of relative coordinates  
#' @export
#' @seealso \code{\link{Profile-class}}, \code{\link{ProfileSet-class}},
#' \code{\link{heightAdjustment}}, \code{\link{calcRelativeCoords}}
#' @examples 
#' # gpsCoordinates = new("GpsCoordinates",
#' #                      address = "../example/gps/p1.txt")
#' data(sinkhole)
#' sinkhole@profiles[[1]]@gpsCoordinates
#' sinkhole@profiles[[1]]@gpsCoordinates@address
#' sinkhole@profiles[[1]]@gpsCoordinates@exact
#' sinkhole@profiles[[1]]@gpsCoordinates@lm
#' sinkhole@profiles[[1]]@gpsCoordinates@relative
#' sinkhole@profiles[[1]]@gpsCoordinates@lmRelative
setClass("GpsCoordinates",
         representation = representation(
           address = "character",
           exact = "data.frame",
           lm = "lm",
           relative = "data.frame",
           lmRelative = "lm"),
         prototype = prototype(
           lm = lm(1~1),
           lmRelative = lm (1~1)))
setMethod("initialize", "GpsCoordinates", 
          function(.Object, address) {
            if(nchar(address) == 0) {
              print("GPS coordinates address is missing.")
            } else {
              .Object@address = address
              
              gpsData <- read.table(file=address, header=T)  
                
              .Object@exact <- data.frame(
                "lat"=gpsData[1],
                "lon"=gpsData[2])
              colnames(.Object@exact) <- c("lat", "lon")
              
              lm <- lm(.Object@exact$lat ~ .Object@exact$lon)
              .Object@lm <- lm
              
              minLat <- min(gpsData[1])
              minLon <- min(gpsData[2])
              
              relativeCoords <- calcRelativeCoords(.Object, minLat, minLon)
              .Object@relative <- relativeCoords
              .Object@lmRelative <- lm(relativeCoords$lat ~ relativeCoords$lon)
            }
            return(.Object)
          })

#' Profile Class
#' 
#' A class to handle a single profile.
#'
#' @slot title title of the profile (e.g. Profile 1).
#' @slot number index of the profile.
#' @slot xyzData object of Xyz Data Class (\code{\link{XyzData-class}}).
#' @slot rawData object of Raw Data Class (\code{\link{RawData-class}}).
#' @slot measurementType type of measurement (e.g. Dipole Dipole, Wenner, ...).
#' @slot gpsCoordinates object of GpsCoordinates Class (\code{\link{GpsCoordinates-class}}).
#' @export
#' @seealso \code{\link{XyzData-class}}, \code{\link{RawData-class}},
#' \code{\link{GpsCoordinates-class}}, \code{\link{plot3dXyz}}
#' @examples 
#' # p1 <- new("Profile",
#' #           title = "Profile 1",
#' #           xyzData = 
#' #             new("XyzData", 
#' #           rawData = 
#' #             new("RawData",
#' #                 address = "../example/rawdata/p1_DipolDipol_SW-NE.dat"),
#' #           measurementType = "DipolDipol",
#' #           gpsCoordinates = 
#' #             new("GpsCoordinates",
#' #                 address = "../example/gps/p1.txt"))
#' #
#' # p1@title
#' # p1@xyzData
#' # p1@rawData
#' # p1@measurementType
#' # p1@gpsCoordinates
#' #
#' # plot3dXyz(p1)
setClass("Profile",
         representation = representation(
           title = "character",
           number = "numeric",
           xyzData = "XyzData",
           rawData = "RawData", 
           measurementType = "character",
           gpsCoordinates = "GpsCoordinates"),
         prototype = prototype(
           number = 0,
           title = "",
           xyzData = NULL,
           rawData = NULL))

#' Profile Set Class
#'
#' A class to handle a collection of many profiles.
#'
#' @slot title title to plot
#' @slot profiles list that contains objects of class Profile (\code{\link{Profile-class}})
#' @slot minLat minimum latitude value of all profiles
#' @slot minLon minimum longitude value of all profiles
#' @slot minData minimum data value of all profiles
#' @slot maxData maximum data value of all profiles
#' @export
#' @seealso \code{\link{Profile-class}}, \code{\link{plot3dXyz}}
#' @examples 
#' # sinkhole <- new("ProfileSet",
#' #                profiles = list(p1, p2, p3),
#' #                title="Sinkhole")
#' 
#' data(sinkhole)
#' plot3dXyz(sinkhole)
setClass("ProfileSet",
         representation = representation(
           profiles = "list",
           title = "character",
           minLat = "numeric",
           minLon = "numeric",
           minData = "numeric",
           maxData = "numeric"))
setMethod("initialize", "ProfileSet",
          function(.Object, profiles=list(), title="",
                   minData=9999999, maxData=0,
                   minLat=100000000000, minLon=100000000000) {
            .Object@profiles <- profiles
            .Object@title <- title         
            
            for (profile in profiles) {
              minDataX <- min(profile@xyzData@seaLevel$val)
              maxDataX <- max(profile@xyzData@seaLevel$val)
              if(minDataX < minData) minData <- minDataX
              if(maxDataX > maxData) maxData <- maxDataX
              
              minLatX <- min(profile@gpsCoordinates@exact$lat)
              minLonX <- min(profile@gpsCoordinates@exact$lon)
              if(minLatX < minLat) minLat <- minLatX
              if(minLonX < minLon) minLon <- minLonX
            }
            
            .Object@minLat <- minLat
            .Object@minLon <- minLon
            .Object@minData <- minData
            .Object@maxData <- maxData
            
            number <- 1
            for(profile in profiles) {
              .Object@profiles[[number]]@number <- number
              
              relativeCoords <- calcRelativeCoords(profile@gpsCoordinates, minLat, minLon)
              .Object@profiles[[number]]@gpsCoordinates@relative <- relativeCoords
              .Object@profiles[[number]]@gpsCoordinates@lmRelative <- lm(relativeCoords$lat ~ relativeCoords$lon)
              
              number <- number + 1
            }
            return(.Object)
          })