### Susumu Tanimura <aruminat@gmail.com>
### Time-stamp: <2015-06-24 19:20:49 umusus>
### Draw very simplified Japan map with prefecture boundaries.

JapanPrefecturesMap <- function(col = NULL, inset = TRUE, ...){
  m <- readShapePoly(system.file("shapes/jpn.shp", package="Nippon")[1],
  proj4string = CRS("+proj=longlat +datum=WGS84"))
  if(inset){
    xy.okinawa <- m@polygons[[47]]@Polygons[[1]]@coords
    xy.okinawa[, 1] <- xy.okinawa[, 1] + 7
    xy.okinawa[, 2] <- xy.okinawa[, 2] + 14
    m@polygons[[47]]@Polygons[[1]]@coords <- xy.okinawa
    labpt <- m@polygons[[47]]@labpt
    labpt[1] <- labpt[1] + 7
    labpt[2] <- labpt[2] + 14
    m@polygons[[47]]@labpt <- labpt
    m@bbox[1, 1] <- 130
    m@bbox[2, 1] <- 31
  }
  plot(m, col = col, ...)
  if(inset){
  lines(x = c(132, 135, 137, 137),
        y = c(38, 38, 40, 43))}
 return(invisible(coordinates(m))) 
}


