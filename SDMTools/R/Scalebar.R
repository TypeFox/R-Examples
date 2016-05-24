#' Scalebar for Projected Maps
#' 
#' \code{Scalebar} adds a distance scalebar onto a projected map. It is not
#' appropriate for geographic projections.
#' 
#' 
#' @param x the x-axis position for the lower left corner of the bar
#' @param y the x-axis position for the lower left corner of the bar
#' @param distance the distance for which the scale bar should represent
#' @param unit the units to report as the scaling
#' @param scale the scaling factor to rescale the distance to a different unit.
#' e.g., if your map is in m and want the scalebar to be in km, use a scale of
#' 0.01
#' @param t.cex the scaling of the font size to be used for the scalebar
#' @return nothing is returned, simply a scalebar is added to a plot.
#' @author Jeremy VanDerWal \email{jjvanderwal@@gmail.com}
#' @examples
#' 
#' 
#' #create a simple object of class 'asc'
#' tasc = as.asc(matrix(1:50,nr=50,nc=50)); print(tasc)
#' 
#' #plot the image
#' image(tasc,axes=FALSE,ann=FALSE)
#' 
#' #add a distance scalebar
#' Scalebar(x=5,y=5,distance=20) #show values in km
#' Scalebar(x=5,y=10,distance=20,unit='m',scale=1000) #show values in meters
#' 
#' 
#' @export 
Scalebar = function(x,y,distance,unit='km',scale=1,t.cex=0.8) {
  xvals = distance*c(0,0.25,0.5,0.75,1)+x
  yvals = c(0,distance/c(30,20,10))+y
  cols <- c("black","white","black","white")
  for (i in 1:4) rect(xvals[i],yvals[1],xvals[i+1],yvals[2],col=cols[i])
  for (i in 1:5) segments(xvals[i],yvals[2],xvals[i],yvals[3])
  labels <- c((xvals[c(1,3)]-xvals[1])*scale,paste((xvals[5]-xvals[1])*scale,unit))
  text(xvals[c(1,3,5)],yvals[4],labels=labels,adj=.5,cex=t.cex)
}
