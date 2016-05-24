#' Draws the contour {f=c}.
#'
#' @param ... An image. Either as a list with components x,y and z or as vectors
#'            x and y and a matrix z of dimensions c(length(x),length(y)).
#' @param level The level of the contour.
#' @param col Color of the contour.
#' @param lty Line type for the contour.
#' @return NULL
#' @export
DrawContour<-function(...,level,col,lty=1){
  C<-contourLines(...,levels=level,nlevels=1)
  if(length(C) == 0) return()
  for(i in 1:length(C))
    lines(C[[i]]$x,C[[i]]$y,pch=20,col=col,lwd=3,lty=lty)
}