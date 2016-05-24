## =============================================================================
## 3-D surfaces
## =============================================================================
# x, y, z, colvar: matrices

surf3D <- function(x, y, z, ..., 
                   colvar = z, phi = 40, theta = 40,
                   col = NULL, NAcol = "white", breaks = NULL,
                   border = NA, facets = TRUE,
                   colkey = NULL, 
                   panel.first = NULL,
                   clim = NULL, clab = NULL, bty = "n",
                   lighting = FALSE, shade = NA, ltheta = -135, lphi = 0,
                   inttype = 1, add = FALSE, plot = TRUE) {

 # check validity, class and dimensionality
  if (! is.matrix(x))
    stop("'x' should be a matrix")
  if (! is.matrix(y))
    stop("'y' should be a matrix")
  if (! is.matrix(z))
    stop("'z' should be a matrix")
  if (ispresent(colvar))
    if (! is.matrix(colvar))
      stop("'colvar' should be a matrix or absent")

  DD <- dim(x)
  if (any (DD != dim(y)) )
    stop("dimension of 'x' not equal to dimension of 'y'")
  if (any (DD != dim(z)) )
    stop("dimension of 'x' not equal to dimension of 'z'")

# check if col or colvar already have the colors to be used

  if (is.character(colvar) & is.matrix(colvar)) {
    col <- colvar
    colvar <- NULL
  }

  if (is.null(col))
    if (is.null(breaks))
      col <- jet.col(100)
    else
      col <- jet.col(length(breaks)-1)

  if (is.null(colvar) & is.matrix(col)) {
    pmat <- persp3Db(x = x, y = y, z = z, col = col, ..., 
             phi = phi, theta = theta, NAcol = NAcol, border = border, 
             facets = facets, panel.first = panel.first,
             bty = bty, lighting = lighting, add = add, plot = plot)
    return(invisible(pmat))
  }

  plist <- initplist(add)

  dot <- splitdotpersp(list(...), bty, lighting, 
    x, y, z, plist = plist, shade, lphi, ltheta, breaks = breaks)
  
  CC <- check.colvar.persp(colvar, z, col, inttype, clim, dot$alpha)
  colvar <- CC$colvar
  col <- CC$col
  
  if (ispresent(colvar)) {

    if (is.null(clim)) 
      clim <- range(colvar, na.rm = TRUE)     
  
    if (dot$clog) {     
      colvar <- log(colvar)
      clim <- log(clim) 
    }

    iscolkey <- is.colkey(colkey, col)
    if (iscolkey) 
      colkey <- check.colkey(colkey)
  
  } else 
    iscolkey <- FALSE

  Extend <- inttype == 2

  if (is.null(plist)) {
    do.call("perspbox", c(alist(x = range(x), y = range(y), 
             z = range(z, na.rm = TRUE),
             phi = phi, theta = theta, plot = plot, 
             colkey = colkey, col = col), dot$persp))
    plist <- getplist()
  }
  breaks <- check.breaks(breaks, col)
  if (is.function(panel.first))
    panel.first(plist$mat)  
           
 # polygons using painters algorithm
  Poly <- paintit(colvar, x, y, z, plist, col, NAcol, clim, border, 
          facets, dot$points$lwd, dot$points$lty, dot, Extend, 
          breaks = breaks)

  if (iscolkey)  
    plist <- plistcolkey(plist, colkey, col, clim, clab, 
          dot$clog, type = "surf3D", breaks = breaks)

  plist <- plot.struct.3D(plist, poly = Poly, plot = plot)  

  setplist(plist)  
  invisible(plist$mat)
}

