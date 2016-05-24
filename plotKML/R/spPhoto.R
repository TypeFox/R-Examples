# Purpose        : Generation of SpatialPhotoOverlay object 
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl);
# Contributions  : ; 
# Dev Status     : Pre-Alpha
# Note           : Combination of the pixmap package and EXIF information [http://www.sno.phy.queensu.ca/~phil/exiftool/TagNames/EXIF.html];


## Generate SpatialPhotoOverlay object:
spPhoto <- function(
   filename,
   obj,
   pixmap,
   exif.info = NULL,
   ImageWidth = 0,
   ImageHeight = 0,
   bands = rep(rep(1, ImageHeight*ImageWidth), 3),
   bbox = c(0,0,3/36000*ImageWidth,3/36000*ImageHeight),
   DateTime = "",
   ExposureTime = "",
   FocalLength = "50 mm", 
   Flash = "No Flash",
   rotation = 0, 
   leftFov = -30,
   rightFov = 30,
   bottomFov = -30,
   topFov = 30,
   near = 50, # in meters;
   shape = c("rectangle", "cylinder", "sphere")[1],
   range = 1000, # m;
   tilt = 90,
   heading = 0,
   roll = 0,
   test.filename = TRUE
   ){

   if(test.filename==TRUE){
     if(!file.exists(filename)){
       if(requireNamespace("RCurl", quietly = TRUE)){
          z <- RCurl::getURI(filename, .opts=RCurl::curlOptions(header=TRUE, nobody=TRUE, transfertext=TRUE, failonerror=FALSE, ssl.verifypeer = FALSE))
          if(!length(x <- grep(z, pattern="404 Not Found"))==0){
            stop(paste("File", filename, "could not be located."))
          } else {
            pixmap <- pixmapRGB(bands, ImageHeight, ImageWidth, bbox = bbox) 
          }
       } else {
         stop('package "RCurl" required but missing')
       }
     } else {
       stop(paste("File", filename, "could not be located."))
     }
   }
    
  # Local copy or in memory
  if(!missing(pixmap)&missing(filename)){
    filename = ""
  }

  ## if missing the coordinate system assume latlon:
  if(!missing(obj)){
    if(is.na(proj4string(obj))) { proj4string(obj) <- CRS(get("ref_CRS", envir = plotKML.opts)) }
  }

  ## if missing EXIF data:
  if(is.null(exif.info)){
    exif.info <- as.list(data.frame(DateTime, ExposureTime, FocalLength, Flash))
  }
  else{ 
    ## try to guess coordinates from EXIF data:
    if(missing(obj)&any(names(exif.info) %in% "GPSLongitude")){
      if(any(names(exif.info) %in% "GPSAltitude")){
        x <- as.numeric(strsplit(exif.info$GPSAltitude, "/")[[1]])
        try(exif.info$GPSAltitude <- ifelse(length(x)>1, x[1]/x[2], x))
      }
      else {
        exif.info$GPSAltitude <- 0
      }
    
    obj <- data.frame(lon=as.numeric(exif.info$GPSLongitude), lat=as.numeric(exif.info$GPSLatitude), alt=as.numeric(exif.info$GPSAltitude))
    coordinates(obj) <- ~lon+lat+alt
    proj4string(obj) <- CRS(get("ref_CRS", envir = plotKML.opts))
    
    }
    else {
      stop("GPS Longitude/Latitude tags not available from the exif.info object.")
    } 
    
    ## correct the ViewVolume:
    exif.info$ImageWidth <- as.numeric(exif.info$ImageWidth)
    exif.info$ImageHeight <- as.numeric(exif.info$ImageHeight)
    asp = exif.info$ImageWidth / exif.info$ImageHeight
    leftFov = leftFov * asp
    rightFov = rightFov * asp
    
    ## format the DateTime field:
    exif.info$DateTime <- format(as.POSIXct(exif.info$DateTime, format="%Y:%m:%d %H:%M:%S", tz="GMT"), "%Y-%m-%dT%H:%M:%SZ")
    
    ## add missing columns:
    if(!any(names(exif.info) %in% "ExposureTime")){
       exif.info$ExposureTime <- ExposureTime
    }
    if(!any(names(exif.info) %in% "FocalLength")){
       exif.info$FocalLength <- FocalLength
    }
    if(!any(names(exif.info) %in% "Flash")){
       exif.info$Flash <- Flash
    }
    
  }
  
  ## Get the heading (if available):
  if(any(names(exif.info) %in% "GPSImgDirection")){
      x <- as.numeric(strsplit(exif.info$GPSImgDirection, "/")[[1]])
      try(exif.info$GPSImgDirection <- ifelse(length(x)>1, x[1]/x[2], x))
      heading = exif.info$GPSImgDirection
  }
  
  ## Photo geometry:
  PhotoOverlay <- as.list(data.frame(rotation, leftFov, rightFov, bottomFov, topFov, near, shape, range, tilt, heading, roll))
  
  ## make a SpatialPhotoOverlay object:
  spPh <- new("SpatialPhotoOverlay", filename = filename, pixmap = pixmap, exif.info = exif.info, PhotoOverlay = PhotoOverlay, sp = obj) 
  return(spPh)    
}

# end of script; 
