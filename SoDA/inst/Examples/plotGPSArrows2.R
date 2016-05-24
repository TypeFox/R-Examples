plotGPSArrows <- function(object, 
                          fraction = 2, head = .5, ...) {
      xy <- geoXY(object@latitude, object@longitude)
      x = xy[,1]; y = xy[,2]
      plot(x, y, asp = 1, type = "n")
      speed <- trackSpeed(cbind(xy,object@elevation), object@time)
      speed <- speed * ( fraction/max(speed, na.rm=TRUE))
      trackArrows(x,  y,  speed, head, ...)
}
