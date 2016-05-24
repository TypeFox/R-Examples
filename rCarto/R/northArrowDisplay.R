northArrowDisplay <-
function(northArrow,northArrowSize){
  if(northArrow==TRUE){
    l<-locator(n=1)
    SpatialPolygonsRescale(layout.north.arrow(2),offset=c(l$x,l$y),
                           scale=northArrowSize,plot.grid=F)
  }
}
