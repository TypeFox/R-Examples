## =============================================================================
## Perspective plot, x, y matrix or vector; z = matrix
## =============================================================================

persp3D <- function(x = seq(0, 1, length.out = nrow(z)), 
                  y = seq(0, 1, length.out = ncol(z)), 
                  z, ..., colvar = z, 
                  phi = 40, theta = 40,
                  col = NULL,  NAcol = "white", breaks = NULL,
                  border = NA, facets = TRUE,
                  colkey = NULL, resfac = 1, 
                  image = FALSE, contour = FALSE, panel.first = NULL,
                  clim = NULL, clab = NULL, bty = "b",
                  lighting = FALSE, shade = NA, ltheta = -135, lphi = 0,
                  inttype = 1, curtain = FALSE, add = FALSE, plot = TRUE){

  plist <- initplist(add)
      
  dot <- splitdotpersp(list(...), bty, lighting, 
    x, y, z, plist = plist, shade, lphi, ltheta, breaks = breaks)
    
 # check dimensionality
  if (! is.vector(x) & length(dim(x)) == 1)
    x <- as.vector(x)
    
  if (! is.vector(y) & length(dim(y)) == 1)
    y <- as.vector(y)

  if (! is.vector(x) & length(dim(x)) > 2)
    stop("'x' should be a vector or a matrix")
    
  if (! is.vector(y) & length(dim(y)) > 2)
    stop("'y' should be a vector or a matrix")

  if ((is.matrix(x) & ! is.matrix(y)) | (is.matrix(y) & ! is.matrix(x)))
    stop ("'x' and 'y' should be of same type, i.e. either a vector or a matrix")

  if (is.matrix(x)) {
    if (any(dim(x) != dim(y)))
      stop ("'x' and 'y' not of same dimension")
    if (curtain)
      stop("cannot combine 'x' a matrix and 'curtain = TRUE'")
  }  
  if (! is.matrix(z)) {
    if (length(z) > 1)
      stop("'z'  should be a matrix or one value") 
    if (is.vector(x))  
      z <- matrix(nrow = length(x), ncol = length(y), data = z)  
    else
      z <- matrix(nrow = nrow(x), ncol = ncol(x), data = z)
  } else if (! is.matrix(x)) {
    if (length(x) != nrow(z))
      stop("'x' should be of length = nrow(z)")
    if (length(y) != ncol(z))
      stop("'y' should be of length = ncol(z)")
  } else
    if (any(dim(z) != dim(x)))
      stop ("'z' and 'x' not of same dimension")

 # resolution    
  if (! is.matrix(x)) {
    if (any(resfac != 1)) {   # change resolution
      res <- changeres(resfac, x, y, z, colvar)
      x <- res$x ; y <- res$y ; z <- res$z
      colvar <- res$colvar
    }
  
 # swap if decreasing
    if (all(diff(x) < 0)) {     
      x <- rev(x)
      if (is.null(dot$persp$xlim)) 
        dot$persp$xlim <- range(x)
      else if (diff(dot$persp$xlim) < 0)
        stop("'persp' expects increasing xlim")  
      if (ispresent(colvar)) 
        colvar <- colvar[nrow(colvar):1, ]
      z <- z[nrow(z):1, ]
    }
    if (all(diff(y) < 0)) {     
      y <- rev(y)
      if (is.null(dot$persp$ylim)) 
        dot$persp$ylim <- range(y)
      else if (diff(dot$persp$ylim) < 0)
        stop("'persp' expects increasing ylim")  
      if (ispresent(colvar)) 
        colvar <- colvar[, (ncol(colvar):1)]
      z <- z[, (ncol(z):1)]
    }
  }

# check if col or colvar already have the colors to be used
  if (is.character(colvar) & is.matrix(colvar)) {
    col <- colvar
    colvar <- NULL
  }

  if (is.null(colvar) & is.matrix(col)) {
    pmat <- persp3Db(x = x, y = y, z = z, col = col, ..., 
             phi = phi, theta = theta, NAcol = NAcol, breaks = breaks,
             border = border, facets = facets, panel.first = panel.first,
             bty = bty, lighting = lighting, shade = shade, ltheta = ltheta,
             lphi = lphi, add = add, plot = plot)
    return(invisible(pmat))
  }
    
  image   <- check.args(image)
  contour <- check.args(contour)
  if (image$add & is.matrix(x))  
      stop("cannot combine 'x' a matrix and 'image'")
  if (contour$add & is.matrix(x))  
      stop("cannot combine 'x' a matrix and 'contour'")

  cv <- colvar

  if (is.null(col) & is.null(breaks))
   col <- jet.col(100)
  else if (is.null(col))
   col <- jet.col(length(breaks)-1)

 # check colvar and colors
  CC <- check.colvar.persp(colvar, z, col, inttype, clim, dot$alpha)
  colvar <- CC$colvar
  col <- CC$col

  Extend <- inttype == 2
  
  if (ispresent(colvar)) {

    if (is.null(clim) & is.null(breaks))
      clim <- range(colvar, na.rm = TRUE)
    else if (is.null(clim))
      clim <- range(breaks, na.rm = TRUE)

    iscolkey <- is.colkey(colkey, col)
    if (iscolkey) 
      colkey <- check.colkey(colkey)
    
    if (dot$clog) {                
      colvar <- log(colvar)
      clim <- log(clim)
    }

  } else 
    iscolkey <- FALSE

  if (is.null(plist)) {
    do.call("perspbox", c(alist(x, y, z,  
                     phi = phi, theta = theta, plot = plot, 
                     colkey = colkey, col = col), dot$persp))
    plist <- getplist()
  }

  breaks <- check.breaks(breaks, col)

  if (is.function(panel.first)) 
    panel.first(plist$mat)         

 # polygon plotting
  if (! is.matrix(x)) { 
    X <- matrix(nrow = nrow(z), ncol = ncol(z), data = x)
    Y <- matrix(nrow = nrow(z), ncol = ncol(z), data = y, byrow = TRUE)
  } else {
    X <- x
    Y <- y
  }  
  lwd <- ifelse (is.null (dot$points$lwd), 1, dot$points$lwd)
  lty <- ifelse (is.null (dot$points$lty), 1, dot$points$lty)

  Poly <- paintit (colvar, X, Y, z, plist, col, NAcol, clim, 
           border, facets, lwd, lty, dot, Extend, breaks = breaks)

  if (curtain) {
    P <- list(x = NULL, y = NULL, col = NULL, border = NULL, 
               lwd = NULL, lty = NULL, proj = NULL, img = list(), 
               isimg = NULL)                    
 
    zmin <- plist$zlim[1]
    Nx <- length(x)
    Ny <- length(y)      
    P <- add.poly(P, 
          cbind(rep(x[1], Ny), rep(x[1], Ny)), cbind(y, y), 
          cbind(rep(zmin, Ny), z[1,]), colvar[1,], 
          col, NAcol, breaks, clim, facets, border, lwd, lty)
    P <- add.poly(P, 
          cbind(rep(x[Nx], Ny), rep(x[Nx], Ny)), cbind(y, y), 
          cbind(rep(zmin, Ny), z[Nx,]), colvar[Nx-1,], 
          col, NAcol, breaks, clim, facets, border, lwd, lty)
    P <- add.poly(P, 
          cbind(x, x), cbind(rep(y[1], Nx), rep(y[1], Nx)), 
          cbind(rep(zmin, Nx), z[,1]), colvar[, 1], 
          col, NAcol, breaks, clim, facets, border, lwd, lty)
    P <- add.poly(P, 
          cbind(x, x), cbind(rep(y[Ny], Nx), rep(y[Ny], Nx)), 
          cbind(rep(zmin, Nx), z[,Ny]), colvar[, Ny-1], 
          col, NAcol, breaks, clim, facets, border, lwd, lty)
    if (! dot$shade$type == "none") {
      P <- color3D(P, plist$scalefac, dot$shade, lighting)
      if (!facets) P$col[] <- "white"
    }
    
   # depth view of the points 
     P$proj   <- project(colMeans(P$x, na.rm = TRUE), 
                         colMeans(P$y, na.rm = TRUE), 
                         colMeans(P$z, na.rm = TRUE), plist)
      
    Poly <- 
      list(x       = cbind(Poly$x, P$x),
           y       = cbind(Poly$y, P$y),               
           z       = cbind(Poly$z, P$z),               
           col     = c(Poly$col, P$col),
           border  = c(Poly$border, P$border),
           proj    = c(Poly$proj, P$proj),
           lwd     = c(Poly$lwd, P$lwd),
           lty     = c(Poly$lty, P$lty),
           isimg   = c(Poly$isimg, P$isimg),
           img     = Poly$img
           )
  class(Poly) <- "poly"
       

  }
   
 # images and contours
  if (image$add) 
    Poly <- XYimage (Poly, image, x, y, z, plist, col, breaks = breaks)

  if (contour$add) 
    segm <- contourfunc(contour, x, y, z, plist, cv = cv, clim = clim)
  else
    segm <- NULL

  if (iscolkey) 
    plist <- plistcolkey(plist, colkey, col, clim, clab, dot$clog, 
      type = "persp3D", breaks = breaks)

  plist <- plot.struct.3D(plist, poly = Poly, segm = segm, plot = plot)  

  setplist(plist)   
  invisible(plist$mat)
}

