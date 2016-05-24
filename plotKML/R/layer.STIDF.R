# Purpose        : Writing of irregular space-time objects to KML
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl);
# Contributions  : Benedikt Graeler (ben.graeler@uni-muenster.de), Pierre Roudier (pierre.roudier@landcare.nz); Edzer Pebesma (edzer.pebesma@uni-muenster.de); 
# Status         : Alpha
# Note           : This method works only with the Space time irregular data frame class objects from the spacetime package;

kml_layer.STIDF <- function(
  obj,
  dtime, 
  ...
  ){

  ## Format the time slot for writing to KML:
  TimeSpan.begin = format(time(obj@time), "%Y-%m-%dT%H:%M:%SZ")
  if(missing(dtime)){  
    TimeSpan.end = format(obj@endTime, "%Y-%m-%dT%H:%M:%SZ")  ## time support BG: all ST* objects now have a slot "endTime"
  } else {
    TimeSpan.end <- format(as.POSIXct(unclass(as.POSIXct(time(obj@time))) + dtime, origin="1970-01-01"), "%Y-%m-%dT%H:%M:%SZ")
  }
  
  # Check the data type:
  if(is(obj@sp, "SpatialPixels")) {
    ## construct stack of rasters:
    r <- brick(lapply(unique(index(obj@time)), 
                      function(x) {
                        raster(obj[,index(obj@time) %in% x, drop=TRUE])
                      }))
    r <- setZ(r, as.character(unique(index(obj@time))))
    if(missing(dtime)){
      dtime <- unique(index(obj@time)) - unique(as.Date(obj@endTime))
      units(dtime) <- "secs"
    }
    kml_layer.RasterBrick(obj = r, dtime=as.numeric(dtime), ...)
  } else {
    if(is(obj@sp, "SpatialPoints")) {
      sp <- SpatialPointsDataFrame(obj@sp, obj@data)
      kml_layer.SpatialPoints(obj = sp, TimeSpan.begin = TimeSpan.begin, TimeSpan.end = TimeSpan.end,  ...)
    } else {
      if(class(obj@sp)=="SpatialPolygons"|class(obj@sp)=="SpatialPolygonsDataFrame"){
        sp <- SpatialPolygonsDataFrame(obj@sp, obj@data)   
        kml_layer.SpatialPolygons(obj = sp, TimeSpan.begin = TimeSpan.begin, TimeSpan.end = TimeSpan.end,  ...)  
      } else {
        if(class(obj@sp)=="SpatialLines"|class(obj@sp)=="SpatialLinesDataFrame"){
          sp <- SpatialLinesDataFrame(obj@sp, obj@data)   
          kml_layer.SpatialLines(obj = sp, TimeSpan.begin = TimeSpan.begin, TimeSpan.end = TimeSpan.end,  ...)
        } else { 
          stop("The 'sp' slot of the ST-object does not extend SpatialPoints*, SpatialLines* or SpatialPolygons*")
        }
      }
    }
  }
}

setMethod("kml_layer", "STIDF", kml_layer.STIDF)
setMethod("kml_layer", "STFDF", function(obj, ...) kml_layer.STIDF(as(obj, "STIDF"), ...))
setMethod("kml_layer", "STSDF", function(obj, ...) kml_layer.STIDF(as(obj, "STIDF"), ...))

# end of script;