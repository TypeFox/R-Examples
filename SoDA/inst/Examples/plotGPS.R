.arrowFrom <- function(u, fraction = .25) {
      n = length(u)
      if(n < 2)
       numeric(0)
      else
      u[-1]*(1-fraction) + u[-n]*fraction
  }

.trackArrows <- function(x, y, z, fraction = .25, col = seq(5)) {
    x0 = .arrowFrom(x, fraction)
    y0 = .arrowFrom(y, fraction)
    x1 = x[-1]
    y1 = y[-1]
    zGroups = cut(z[!is.na(z)], length(col))
    ## compute the average line length
    delta = sqrt(mean((x1-x0)^2 + (y1-y0)^2, na.rm = TRUE))
    ## and convert it to inches for arrows()
    ## Assumes that plot coordinates were set up with asp=1 or the
    ## equivalent
    delta = delta * (par("pin")[1]/diff(range(x, na.rm = TRUE)))
    arrows(x0, y0, x1, y1, delta, col = zGroups)
}
    
    
plotGPSArrows <- function(latitude, longitude, elevation,
                          fraction = .25, add = FALSE, ...) {
      xy <- geoXY(latitude, longitude)
      x = xy[,1]; y = xy[,2]
     ## set up plot
      if(!add)
          plot(range(x, na.rm = TRUE), range(y, na.rm = TRUE),
                  asp = 1, type = "n", ...)
      .trackArrows(x,  y, elevation, fraction)
}
