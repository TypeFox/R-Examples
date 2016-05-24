## ###########################################################
## PURPOSE: Put color bar legend on the right side of an existing
##    plot.
##
## INPUT:
##    zlim: Gives the range of z values to which the colors specified
##       in 'col' are assigned.
##    col: Gives the range of colors to use.  To keep multiple
##       plots consistent in terms of the colors assigned to various
##       values, keep 'zlim' and 'col' the same for each of the plots
##       and the legend.
##
## RESULT: Puts vertical color bar legend to the right of the plot.
##
## ASSUMES:
##    1. Originally, this function assumed that the user had already
##       finished making the main portion of the plot; i.e., it
##       required that the user should ADD THE LEGEND LAST.  While the
##       author believes this issue has been fixed, the author still
##       recommends that the user avoid changing the main plot after
##       using this legend command.  This function alters the par()
##       settings to draw the legend, and trying to reset them
##       properly is not straightforward.
##    2. This function works best when there is only one plot on the
##       device.
##
## 
## NOTES:  As noted in the help for image.plot() in the "fields"
##    package, putting a legend on a plot is harder than it might seem.
##    The user may have to experiment with this function a bit to get
##    it to work well for a specific application.  The user may also
##    want to try the previously mentioned image.plot() function in
##    the fields package (to just add the legend, use "legend.only=T").
##
##
## REVISION HISTORY:
##   Prototype: Jenise Swall, 2010-04-09
##
##   2010-06-16 (JLS): Altered function so that, when exiting, we do
##      not try to reset the par settings to what they were when we
##      entered the function.  This was causing sporadic errors that
##      are as yet unexplained (at least to the writer of this function).
##   2011-03-01 (JLS): Altered/clarifed comments.
##   2012-05-20 (JLS): Altered function so that, when exiting, we
##      reset par settings so that changes can still be made to the
##      main plot (not the legend).  These changes were closely based
##      on code from image.plot() in the "fields" package.
##   2013-01-11 (JLS): Fixed bug in determining the cut/endpoints to
##      go with each color.  Previously, the zlim bounds were shown as
##      though they were the midpoints of the range for their color,
##      rather than the endpoints.
## ###########################################################
vertical.image.legend <- function(zlim, col){


  ## Get the current par information.  We restore these settings
  ## before we leave this function.
  starting.par.settings <- par(no.readonly=TRUE)


  ## Find information about the overall size of the figure and the
  ## margins included with that figure.
  mai <- par("mai")
  fin <- par("fin")
  

  ## Find total amount of figure space left over after the main plot
  ## is drawn. This really boils down to the area to the right of the
  ## plot, but under the top margin (which may hold a title) and above
  ## the bottom margin (which may hold a label for the x-axis).
  x.legend.fig <- c( 1.0-(mai[4]/fin[1]), 1.0 )
  y.legend.fig <- c( mai[1]/fin[2], 1.0-(mai[3]/fin[2]) )
  
  ## Now, give the portion of this area which can be used for the
  ## actual legend strip.  This means we need to leave a litle space
  ## between the strip and the main plot, and a bigger space for the
  ## labels to the right of the strip.  In the y direction, we can use
  ## all the space available.
  x.legend.plt <- c( x.legend.fig[1]+(0.08*(x.legend.fig[2]-x.legend.fig[1])),
                    x.legend.fig[2]-(0.6*(x.legend.fig[2]-x.legend.fig[1])) )
  y.legend.plt <- y.legend.fig
  

  ## Find cut/endpoints.  The lower zlim is the bottom endpoint; the
  ## upper zlim element is the lower endpoint.  Then, there are
  ## length(col)-1 cutpoints between these to get length(col)
  ## intervals.
  cut.pts <- seq(zlim[1], zlim[2], length=length(col)+1)
  ## Find the midpoint for each interval.  These are the values that
  ## will cause each color to be plotted in the color bar.
  z <- ( cut.pts[1:length(col)] + cut.pts[2:(length(col)+1)] ) / 2
  

  par( new=TRUE, pty="m", plt=c(x.legend.plt, y.legend.plt) )
  ## par( new=T, xpd=T, pty="m", plt=c(x.legend.plt, y.legend.plt) )
  image(x=1, y=z, z=matrix(z, nrow=1, ncol=length(col)),
        col=col, xlab="", ylab="", xaxt="n", yaxt="n")
  axis(4, mgp = c(3, 0.2, 0), las = 2, cex.axis=0.5, tcl=-0.1)
  box()


  ## Return par settings to what they were when we entered this
  ## function.  This is closely based on what is done in at the end
  ## of the image.plot() function in the "fields" package.
  mfg.settings <- par()$mfg
  par(starting.par.settings)
  par(mfg=mfg.settings, new=FALSE)
}
