##' Returns device-coordinates and plot-region
##'
##' @title Returns device-coordinates and plot-region
##' @return A \code{list} with elements
##'  \item{dev.x1}{Device: Left x-coordinate}
##'  \item{dev.x2}{Device: Right x-coordinate}
##'  \item{dev.y1}{Device Bottom y-coordinate}
##'  \item{dev.y2}{Device Top y-coordinate}
##'  \item{fig.x1}{Plot: Left x-coordinate}
##'  \item{fig.x2}{Plot: Right x-coordinate}
##'  \item{fig.y1}{Plot: Bottom y-coordinate}
##'  \item{fig.y2}{Plot: Top y-coordinate}
##' @author Klaus K. Holst
##' @export
##' @keywords hplot
`devcoords` <-
function() {
  cc <- par("usr") ## extremes of coordinates of plotting region (x1,x2,y1,y2)
  plotinch <- par("pin") ## Plot dimensions (width,height) in inches
  margininch <- par("mai") ## Margin sizes in inches (bottom, left, top ,right)
  plotlenX <- cc[2]-cc[1]
  unitinchX <- plotlenX/plotinch[1]
  plotlenY <- cc[4]-cc[3]
  unitinchY <- plotlenY/plotinch[2]
  deviceXleft <- cc[1]-unitinchX*margininch[2]
  deviceXright <- cc[2]+unitinchX*margininch[4]
  deviceYtop <- cc[4]+unitinchY*margininch[3]
  deviceYbottom <- cc[3]-unitinchY*margininch[1]
  return(list(dev.x1=deviceXleft, dev.x2=deviceXright, dev.y1=deviceYbottom, dev.y2=deviceYtop, fig.x1=cc[1], fig.x2=cc[2], fig.y1=cc[3], fig.y2=cc[4]))
}
