plotGPSArrows <- function(object, 
                          fraction = .75, head = .5, ...) {
      xy <- geoXY(object@latitude, object@longitude)
      x = xy[,1]; y = xy[,2]
      plot(x, y, asp = 1, type = "n")
      trackArrows(x,  y, fraction, head, ...)
}
