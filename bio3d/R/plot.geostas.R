plot.geostas <- function(x, at=seq(0, 1, 0.1), main="AMSM with Domain Assignment",
                         col.regions=rev(heat.colors(200)), 
                         margin.segments=x$grps, ...) {
  
  plot.dccm(x$amsm, at=at, main=main, 
            col.regions=col.regions,
            margin.segments=margin.segments, ...)
  
}
