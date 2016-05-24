## #########################################################
## PURPOSE: Given a list of points' coordinates and the values
##    observed at those points, return a scatterplot with points
##    located as specified by the coordinates and coded by color
##    and/or size to represnt the observed value at the location.
##    This code is basically a wrapper for a call to the function
##    points.geodata() in the geoR package.
##
##
## INPUTS:
##    x, y: Coordinates of locations at which observations were made.
##    z: Values observed at the locations whose coordinates are
##       given by (x, y).
##    zlim: Minimum and maximum value of 'z' to which to assign the
##        two most "extreme" colors in the 'col' argument.  Default
##        is to use the range of z.  This is very much like the 'zlim'
##        argument to the image() function.
##    add: If FALSE, put points onto a new plot (default). If TRUE,
##       adds points to a pre-existing plot.        
##    col: Color range to use for the points, with the first color
##        assigned to zlim[1] and last color assigned to zlim[2].
##    xlab, ylab: The labels to assign to the x- and y-axes.  If not
##       specified by the user, these default to the expressions the
##       user named as parameters x and y.  This behavior is similar
##       to that for plot().
##    pch: The point symbol to use.  Possible values are 21, 22, 23,
##         24, 25.  This is because points.geodata() requires these
##         points, which have outlines around them.  Default is a
##         circle (pch=21).
##    cex.min, cex.max: These control the minimum and maximum amounts
##       to shrink/expland the points, based on the value of z.  By
##       default, these are both set to one, which makes all the points
##       the same size.  For more information, see the help page for
##       points.geodata().
##    symbol.border.col: This controls the color of the border around
##       the plotting symbol.  By default, it is black.  If a border is
##       not desired, use symbol.border.col="transparent".
##    plt.beyond.zlim: If TRUE (and if zlim is specified by the user),
##       z values beyond the limits given in zlim are plotted.
##       Values below zlim[1] are plotted in the same color as zlim[1];
##       values above zlim[2] are plotted in the same color as zlim[2].
##    ...: Any other parameters the user adds will be passed to the
##        plot() function if add=F, and may include options for axis
##        labels, titles, etc.
##
##
## RETURNS:  A scatterplot with points at (x, y).  These points are
##    colored according to the correspoinding value of z and the colors
##    specified in col.  They are sized according to the corresponding
##    value of z and the minimum and maximum sizes specified by 'cex.min'
##    and 'cex.max'.
##
## ASSUMES:
##   Availability of package geoR for the key function points.geodata().
##
##
## REVISION HISTORY:
##   Prototype: Jenise Swall, 2010-04-09
##
##   2011-02-28 (JLS): Changed the name of this function from
##      draw.color.scatterplot.w.geoR to plot3d.points.
##      Adjusted/clarified comments.
##
##   2011-10-21 (JLS): When zlim is provided by the user and
##      plt.beyond.zlim=TRUE, allow plotting of points which have z
##      values outside the bounds given by zlim.  For z values
##      exceeding zlim[2], use the color corresponding to zlim[2]; for
##      z values below zlim[1], use the color for zlim[1].   Also,
##      added explicit parameters xlab and ylab.  This makes it easier
##      to test whether the user has provided xlab and ylab, and, if
##      not, allows these labels to default to the expressions the
##      user named as parameters x and y.  This is similar to the
##      default behavior for plot().
## #########################################################
plot3d.points <- function(x, y, z, zlim=range(z, na.rm=TRUE), add=FALSE,
                          col=heat.colors(12), xlab, ylab, pch=21,
                          cex.min=1, cex.max=1,
                          symbol.border.col="black",
                          plt.beyond.zlim=FALSE, ...){

  ## Pch values must be 21, 22, 23, 24, or 25, as required by
  ## points.geodata().
  if (!(pch %in% c(21, 22, 23, 24, 25)))
    stop("Value of pch must be 21, 22, 23, 24, or 25")


  ## If plt.beyond.zlim=TRUE and zlim is given, we set z values
  ## which are greater than zlim[2] to equal zlim[2] and z values
  ## which are less than zlim[1] to equal zlim[1].  This ensures that
  ## these values are plotted, when they would otherwise be blank.  If
  ## plt.beyond.zlim=TRUE and zlim is missing, then zlim will be
  ## determined according to the total range of the z's and
  ## plt.beyond.zlim can be set to FALSE.
  if (plt.beyond.zlim){
    if (missing(zlim)){
      warning("plt.beyond.zlim=TRUE is not a valid option if zlim argument is not provided, changing to plt.beyond.zlim=FALSE")
      plt.beyond.zlim <- FALSE
    }
    else{
      ## Values lower than zlim[1] set to zlim[1].
      z[z < zlim[1]] <- zlim[1]
      ## Values higher than zlim[2] set to zlim[2].
      z[z > zlim[2]] <- zlim[2]
    }
  }
    

  ## Divide the total range up into as many bins as you have colors.
  ## The object formed below contains the endpoints for each bin, so
  ## it has to have length of 1 plus the number of bins.
  how.to.divide <- seq(zlim[1], zlim[2], length=1+length(col))


  ## Note that points.geodata() will set the plot axes so that they
  ## have the same length.  While this makes the perspective better,
  ## it also can result in a lot of white space in the plot.  So, we
  ## set up the plotting region first (using the R default), and then
  ## add the points with the points.geodata function.  This has the
  ## added advantage of letting us control the symbol border color,
  ## which cannot be done in points.geodata unless you're adding the
  ## points to a plot. 
  if (!add){
    ## If the user doesn't specify xlab and ylab, assign to them
    ## the expressions the user named as parameters x and y.  This
    ## is similar to the default behavior for plot().
    if (missing("xlab"))
      xlab <- deparse(substitute(x))
    if (missing("ylab"))
      ylab <- deparse(substitute(y))
    plot(x, y, type="n", xlab=xlab, ylab=ylab, ...)
  }
  ## Make x, y, z into a geoR geodata object.
  points.geodata(coords=cbind(x,y), data=z, col.seq=col,
                 pt.divide=how.to.divide, cex.min=cex.min, cex.max=cex.max,
                 pch.seq=pch, add.to.plot=TRUE, col=symbol.border.col)
}
