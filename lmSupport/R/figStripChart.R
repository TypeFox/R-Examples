figStripChart <- function(x, side=1, sshift=0.3, col='gray', pch=15, cex= 0.2, adjoffset=1) {
# coded by Dr Steven J. Murdoch (http://www.cl.cam.ac.uk/~sjm217/index.html#contact)
# Draw a stripchart on an axis, showing marginal frequency
# TODO: Does not handle log axes well 
# x:    the data from which the plots are to be produced.
# side: as in axis(). default is x-axis
# col of points
# pch of points
# cex of points

 if (side==1){
   flip <- 1# Does the outside of the plot have larger or smaller vales
   yaxis <- FALSE # Are we on the yaxis
   parside <-3  # Relevant index of par("usr")
 }
 else if (side==2) {
   flip <- 1
   yaxis <- TRUE
   parside <-1
 }
 else if (side==3) {
   flip <- -1
   yaxis <- FALSE
   parside <-4
 }
 else if (side==4) {
   flip <- -1
   yaxis <- TRUE
   parside <-2
 }

 # Axis position
 base<-par("usr")[parside]

 # Width and height in user units
 plotwidth <- diff(par("usr")[1:2])
 plotheight <- diff(par("usr")[3:4])
 
 shift <- par("pin")[1]*0.003*flip  # Shift for the q2 and q3 axis from the base (in inches)
 gap <- par("pin")[1]*0.003  # Gap for the median
 meanshift <- par("cin")[1]*0.5*flip  # Shift for the mean pointer away from the axis
 stripshift <- par("cin")[1]*sshift*flip  # Shift away from the q2 and q3 axis for the stripchart

 # Scale lengths so both axes are equal on output device
 if (yaxis) {
   shift <- shift/par("pin")[1]*plotwidth
   meanshift <- meanshift/par("pin")[1]*plotwidth
   stripshift <- stripshift/par("pin")[1]*plotwidth
   gap <- gap/par("pin")[2]*plotheight
 } else {
   shift <- shift/par("pin")[2]*plotheight
   meanshift <- meanshift/par("pin")[2]*plotheight
   stripshift <- stripshift/par("pin")[2]*plotheight
   gap <- gap/par("pin")[1]*plotwidth
 }

 # If vertical, stripchart assumes offset is a factor of character
 # width, if horizontal, character height (bug?). So correct for this
 if (yaxis)
   offset=flip*par("cin")[2]/par("cin")[1] * adjoffset
 else
   offset=flip* adjoffset

 # Don't clip the chart
 oldxpd <- par(xpd = TRUE)
 on.exit(par(oldxpd))

 stripchart(x, method="stack", vertical=yaxis, offset=offset, pch=pch,
            cex=cex, add=TRUE, at=base+shift+stripshift, col=col)
}