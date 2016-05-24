# ------------------------------------------------------------------------------
# Internal function 'surface.equal'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------
surface.equal <- function(x, data, nrow, ncol, verbose) {

  if (verbose){
    begTime <- Sys.time(); fn <- match.call()[[1]]
    message(fn, ": re-distribute the population data ...")
  }
  
  if (!inherits(x, "SpatialPolygons"))     
    stop(paste("'x' must be a SpatialPolygons object to",
               "use the \"equal\" smoothing option"), call. = FALSE)

  xy <- bbox(x)
  xmn <- xy[1,1]; xmx <- xy[1,2]
  ymn <- xy[2,1]; ymx <- xy[2,2]
  
  # Create a nrow * ncol grid point data set
  xx <- seq(xmn, xmx, length.out = ncol)
  yy <- seq(ymn, ymx, length.out = nrow)
  coords <- expand.grid(xx, yy)
  colnames(coords) <- c("x", "y")
  
  # Points belong to which polygon?  
  spcoords <- SpatialPoints(coords, proj4string = x@proj4string)
  if (is(x, "SpatialPolygonsDataFrame"))
    x <- as(x, "SpatialPolygons")
  polygonIDs <- sp::over(spcoords, x)

  # Remove points that are outside of the polygons
  outside <- which(is.na(polygonIDs))
  if (length(outside) > 0) {
    coords <- coords[-outside,]
    polygonIDs <- polygonIDs[-outside]  
  }
    
  # Now re-distribute the population within each polygon
  cellPerPolygon <- as.vector(table(polygonIDs))
  values <- matrix(NA, nrow = nrow(coords), ncol = ncol(data))
  for (i in 1:ncol(data)) {
    counts <- data[polygonIDs, i]
    counts <- counts / cellPerPolygon[polygonIDs]
    values[,i] <- counts
  }

  if (verbose){
    tt <- as.numeric(difftime(Sys.time(), begTime, units = "sec"))
    message(fn, ": done! [", tt, " seconds]")
  }

  colnames(coords) <- c("x", "y")
  colnames(values) <- colnames(data)
  list(coords = coords, data = values, id = polygonIDs)
}
