## ###########################################################
## PURPOSE: This code produces a ragged image plot.  This is in
##   contrast to the standard image() function, which assumes that
##   there is a known response value 'z' for every combination of
##   the elements of'x' and 'y', i.e. that there is a complete
##   rectangular grid, or image.  A ragged image plot is a variant
##   of the regular image plot which is not complete across the
##   entire rectangle.  The user specifies vectors 'x', 'y', and 'z',
##   such that 'x' and 'y' identify a portion of the grid. This 
##   function maps the vector 'z' vector onto a matrix of the type
##   needed for the image() function, but has NA entries for
##   combinations of x and y that are not listed.  The NA values are
##   not plotted by image(), so a ragged image will appear.
##
##
## INPUT:
##    x, y: Coordinates of grid cell centers for which response values
##          'z' are available.
##    z: Response values recorded at the grid cell centers whose
##       coordinates are given by ('x', 'y').
##    zlim: Minimum and maximum value of 'z' to which to assign the
##        two most "extreme" colors in the 'col' argument (col[1] and
##        col[length(col)].  Default is to use the range of z.  This
##        is very much like the 'zlim' argument to the image function.
##    add: If FALSE, the ragged image will begin a new plot. If TRUE,
##       adds ragged image to a pre-existing plot.        
##    col: Color range to use for the ragged image, with the first
##       color assigned to zlim[1] and last color assigned to zlim[2].
##    xlab, ylab: The labels to assign to the x- and y-axes.  If not
##       specified by the user, these default to the expressions the
##       user named as parameters x and y.  This behavior is similar
##       to that for image().
##    plt.beyond.zlim: If TRUE (and if zlim is specified by the user),
##       z values beyond the limits given in zlim are plotted.
##       Values below zlim[1] are plotted in the same color as zlim[1];
##       values above zlim[2] are plotted in the same color as zlim[2].
##    ...: Any additional parameters to be passed to image function,
##         if add=FALSE.
##
##
## RETURNS: A ragged image (i.e., a portion of an image).
##
##
## NOTES: This function is slow if 'x', 'y', and 'z' are long vectors.
##
##
## REVISION HISTORY:
##   Prototype: Jenise Swall, 2011-03-01 (based on Jenise's function
##      image.irreg.xgrid).
##
##   2011-10-21 (JLS): When zlim is provided by the user and
##      plt.beyond.zlim=TRUE, allow plotting of points which have z
##      values outside the bounds given by zlim.  For z values
##      exceeding zlim[2], use the color corresponding to zlim[2]; for
##      z values below zlim[1], use the color for zlim[1].  Also,
##      added explicit parameters xlab and ylab.  This makes it easier
##      to test whether the user has provided xlab and ylab, and, if
##      not, allows these labels to default to the expressions the
##      user named as parameters x and y.  This is similar to the
##      default behavior for image().
## ###########################################################
ragged.image <- function(x, y, z, zlim=range(z, na.rm=TRUE), add=FALSE,
                         col=heat.colors(12), xlab, ylab,
                         plt.beyond.zlim=FALSE, ...){

  ## If any elements of x or y (location information) are missing,
  ## then exit the function.
  if ( sum(is.na(x)) > 0 )
    stop("No missing values allowed in 'x'.")
  if ( sum(is.na(y)) > 0 )
    stop("No missing values allowed in 'y'.")
  

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
    
  
  ## Construct a regular grid, making sure that it includes the
  ## complete range of x and y values in the original data.
  x.uniq <- sort(unique(x))
  y.uniq <- sort(unique(y))
  grid.coords <- expand.grid(x.uniq, y.uniq)


  ## Create a new vector, which has the original data in right places
  ## to correspond to the regular grid, and has NAs in the regular
  ## grid spaces for which we don't have data.
  z.surf.vec <- rep(NA, nrow(grid.coords))
  match.locs <- match(paste(grid.coords[,1], grid.coords[,2], sep="-"),
                      paste(x, y, sep="-"))
  z.surf.vec <- z[match.locs]
  z.surf.mat <- matrix(z.surf.vec, ncol=length(y.uniq))
  rm(match.locs, z.surf.vec)

  
  ## If add=TRUE, add the image to the existing plot.
  if (!add){

    ## If the user doesn't specify xlab and ylab, assign to them
    ## the expressions the user named as parameters x and y.  This
    ## is similar to the default behavior for image().
    if (missing("xlab"))
      xlab <- deparse(substitute(x))
    if (missing("ylab"))
      ylab <- deparse(substitute(y))
    image(x.uniq, y.uniq, z.surf.mat, zlim=zlim, col=col,
          xlab=xlab, ylab=ylab, ...)
  }
  ## Otherwise, this image begins a new plot.
  else
    image(x.uniq, y.uniq, z.surf.mat, zlim=zlim, col=col, add=TRUE)
}
