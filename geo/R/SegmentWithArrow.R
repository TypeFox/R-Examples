#' Plot line segment with arrow at the end.
#' 
#' Plot line segment with arrow at the end.
#' 
#' 
#' @param pos Segment positions.
#' @param angle (Half-)Angle (< 90) determining arrow sharpness.
#' @param size Size of arrow.
#' @param col Color of arrow.
#' @param lwd Line width of segment.
#' @note Needs further checking and elaboration.
#' @keywords aplot
#' @export SegmentWithArrow
SegmentWithArrow <-
function(pos,angle=15,size=0.2,col="blue",lwd=2){
  geopar <- getOption("geopar")
  plt.size <- geopar$gpar$pin
  dist <- arcdist(pos$lat[1],pos$lon[1],pos$lat[2],pos$lon[2],scale="nmi")
  arrowsize <- diff(geopar$origin$lat)*size/geopar$gpar$pin[2]*60
  rat <- min(c(arrowsize/dist,0.5))
  tmp <- pos
  tmp[1,] <- tmp[2,]+rat*(tmp[1,]-tmp[2,])
  tmp <- Proj(tmp)
  dx <- diff(tmp$x)
  dy <- diff(tmp$y) 
  rat <- tan(angle*pi/180)
  tmp1 <- data.frame(x=tmp$x[c(1,2,1)],y=tmp$y[c(1,2,1)])
  tmp1$y[1] <- tmp1$y[1]-rat*dx
  tmp1$x[1] <- tmp1$x[1]+rat*dy
  tmp1$y[3] <- tmp1$y[3]+rat*dx
  tmp1$x[3] <- tmp1$x[3]-rat*dy
  tmp1 <- invProj(tmp1)
  if(lwd > 0) geolines(pos,col=col,lwd=lwd)
  geopolygon(tmp1,col=col) 
}

