## =============================================================================
## Perspective plot, x, y matrix or vector; z = matrix
## =============================================================================

ribbon3D <- function(x = seq(0, 1, length.out = nrow(z)), 
                   y = seq(0, 1, length.out = ncol(z)), 
                   z, ..., 
                   colvar = z, phi = 40, theta = 40,
                   col = NULL,  NAcol = "white", breaks = NULL,
                   border = NA, facets = TRUE,
                   colkey = NULL, resfac = 1, 
                   image = FALSE, contour = FALSE, panel.first = NULL,
                   clim = NULL, clab = NULL, bty = "b", 
                   lighting = FALSE, shade = NA, ltheta = -135, lphi = 0,
                   space = 0.4, along = "x", 
                   curtain = FALSE, add = FALSE, plot = TRUE) {

  plist <- initplist(add)

  if (any(space > 0.9))
    stop("'space' too large, should be smaller than or equal to 0.9")
  else if (any(space < 0.1))
    stop("'space' cannot be smaller than 0.1")
  space <- rep(space, length.out = 2) # in x- and y
  
 # input check
  if (length(grep("x", along)) == 0 &  length(grep("y", along)) == 0)
    stop ("'along' should contain at least one of 'x' or 'y'")

  if (! is.vector(x) & length(dim(x)) != 1)
    stop("'x' should be a vector")
    
  if (! is.vector(y) & length(dim(y)) != 1)
    stop("'y' should be a vector")

  if (length(x) != nrow(z))
    stop("'x' should be of length = nrow(z)")

  if (length(y) != ncol(z))
    stop("'y' should be of length = ncol(z)")

  if (any(resfac != 1)) {   # change resolution
    res <- changeres(resfac, x, y, z, colvar)
    x <- res$x ; y <- res$y ; z <- res$z
    colvar <- res$colvar
  }

  rx <- range(x)
  ry <- range(y)

  if (length(grep("x", along)) > 0) {
    dY <- 0.5*(1 - space[2]) * diff(y)
    dY <- c(dY[1], dY, dY[length(dY)])
    ry <- ry + c(- dY[1], dY[length(dY)])
  }
   
  if (length(grep("y", along)) > 0) {
    dX <- 0.5*(1 - space[1]) * diff(x)
    dX <- c(dX[1], dX, dX[length(dX)])
    rx <- rx + c(- dX[1], dX[length(dX)])
  }

  dot <- splitdotpersp(list(...), bty, lighting,
    rx, ry, z, plist = plist, shade, lphi, ltheta, breaks = breaks)

 # swap if decreasing
  if (! is.matrix(x) & all(diff(x) < 0)) {    
    if (is.null(dot$persp$xlim)) 
      dot$persp$xlim <- rev(range(x))
    x <- rev(x)
    if (ispresent(colvar)) 
      colvar <- colvar[nrow(colvar):1, ]
    z <- z[nrow(z):1, ]
  }
 
  if (! is.matrix(y) & all(diff(y) < 0)) {    
    if (is.null(dot$persp$ylim)) 
      dot$persp$ylim <- rev(range(y))
    y <- rev(y)
    if (ispresent(colvar)) 
      colvar <- colvar[, (ncol(colvar):1)]
    z <- z[, (ncol(z):1)]
  }
  image   <- check.args(image)
  contour <- check.args(contour)
  if (contour$add) 
    cv <- colvar

  if (is.null(col) & is.null(breaks))
   col <- jet.col(100)
  else if (is.null(col))
   col <- jet.col(length(breaks)-1)

  breaks <- check.breaks(breaks,col)
  CC <- check.colvar.2(colvar, z, col, clim, dot$alpha)
  colvar <- CC$colvar
  col <- CC$col

  if (ispresent(colvar)) {

    if (is.null(clim)) 
      clim <- range(colvar, na.rm = TRUE)
  
    if (dot$clog) {                    
      colvar <- log(colvar)
      clim   <- log(clim)
    }
  
    iscolkey <- is.colkey(colkey, col)     
    if (iscolkey) 
      colkey <- check.colkey(colkey)
  } else
    iscolkey <- FALSE
  rx <- range(x)
  ry <- range(y) 
  
  if (is.null(plist)) {
    do.call("perspbox", c(alist(x = range(x), y = range(y), 
             z = range(z, na.rm = TRUE),
             phi = phi, theta = theta, plot = plot, colkey = colkey, col = col), 
             dot$persp))
    plist <- getplist()
  }

  if (is.function(panel.first)) 
    panel.first(plist$mat)         

 # draw ribbons as polygons 
  shade <- dot$shade$shade
  if (is.null(dot$shade$shade))
    dot$shade$shade <- NA

  Nx <- dim(z) [1]
  Ny <- dim(z) [2]

  lwd <- dot$points$lwd
  if (is.null(lwd)) 
    lwd <- 1
  lty <- dot$points$lty
  if (is.null(lty)) 
    lty <- 1

  Poly <- list(x = NULL, y = NULL, col = NULL, border = NULL, 
               lwd = NULL, lty = NULL,  
               img = NULL, alpha = NULL, proj = NULL)                    

  zmin <- min(plist$zlim[1], min(z, na.rm = TRUE))

  if (curtain & zmin == min(z, na.rm = TRUE))
    zmin <- as.double(zmin - diff(range(plist$zlim)) * 1e-6)  #  to avoid triangle rather than quad

  if (length(grep("x", along)) > 0) {
    X <- cbind(x, x)
    
    for (i in 1 : Ny) {
      ind <- length(Poly$col) + 1
      Y  <- cbind(rep(y[i]+dY[i+1], Nx), rep(y[i]-dY[i], Nx))  
      Poly <- add.poly(Poly, X, Y,  cbind(z[,i], z[,i]), 
                       colvar[,i], col, NAcol, breaks,
                       clim, facets, border)
      if (curtain) {
        Poly <- add.poly(Poly, X, 
          cbind(rep(y[i]-dY[i], Nx), rep(y[i]-dY[i], Nx)),  
          cbind(rep(zmin, Nx), z[,i]), colvar[,i], 
          col, NAcol, breaks, clim, facets, border)

        Poly <- add.poly(Poly, X, 
          cbind(rep(y[i]+dY[i+1], Nx), rep(y[i]+dY[i+1], Nx)),  
          cbind(rep(zmin, Nx), z[,i]), colvar[,i], 
          col, NAcol, breaks, clim, facets, border)

        ye1 <- y[i]-dY[i]
        ye2 <- y[i]+dY[i+1]

        Poly <- 
          list(x      = cbind(Poly$x, c(x[1], x[1], x[1], x[1], NA)),
               y      = cbind(Poly$y, c(ye1, ye2, ye2, ye1, NA)),               
               z      = cbind(Poly$z, c(zmin, zmin, z[1,i], z[1,i], NA)),               
               col    = c(Poly$col, Poly$col[ind]),
               border = c(Poly$border, Poly$border[ind]),
               img    = Poly$img)

        ind <- length(Poly$col)
        
        Poly <- 
          list(x      = cbind(Poly$x, c(x[Nx], x[Nx], x[Nx], x[Nx], NA)),
               y      = cbind(Poly$y, c(ye1, ye2, ye2, ye1, NA)),               
               z      = cbind(Poly$z, c(zmin, zmin, z[Nx,i], z[Nx,i], NA)),               
               col    = c(Poly$col, Poly$col[ind]),
               border = c(Poly$border, Poly$border[ind]),
               img    = Poly$img)
      }
    }
  }
 
  if (length(grep("y", along)) > 0) {
    Y <- cbind(y, y)

    for (i in 1 : Nx) {
      ind <- length(Poly$col) + 1
      X  <- cbind(rep(x[i]+dX[i+1], Ny), rep(x[i]-dX[i], Ny))  
      Poly <- add.poly(Poly, X, Y,  cbind(z[i,], z[i,]), 
                       colvar[i,], col, NAcol, breaks,
                       clim, facets, border)
      if (curtain) { 
        Poly <- add.poly(Poly, 
          cbind(rep(x[i]-dX[i], Ny), rep(x[i]-dX[i], Ny)), Y, 
          cbind(rep(zmin, Ny), z[i,]), colvar[i,], 
          col, NAcol, breaks, clim, facets, border)

        Poly <- add.poly(Poly,  
          cbind(rep(x[i]+dX[i+1], Ny), rep(x[i]+dX[i+1], Ny)), Y,
          cbind(rep(zmin, Ny), z[i,]), colvar[i,], 
          col, NAcol, breaks, clim, facets, border)
      
        xe1 <- x[i]-dX[i]
        xe2 <- x[i]+dX[i+1]

        Poly <- 
          list(x      = cbind(Poly$x, c(xe1, xe2, xe2, xe1, NA)),               
               y      = cbind(Poly$y, c(y[1], y[1], y[1], y[1], NA)),
               z      = cbind(Poly$z, c(zmin, zmin, z[i,1], z[i,1], NA)),               
               col    = c(Poly$col, Poly$col[ind]),
               border = c(Poly$border, Poly$border[ind]),
               img    = Poly$img)
      
        ind <- length(Poly$col)
        
        Poly <- 
          list(x      = cbind(Poly$x, c(xe1, xe2, xe2, xe1, NA)),               
               y      = cbind(Poly$y, c(y[Ny], y[Ny], y[Ny], y[Ny], NA)),
               z      = cbind(Poly$z, c(zmin, zmin, z[i,Ny], z[i,Ny], NA)),               
               col    = c(Poly$col, Poly$col[ind]),
               border = c(Poly$border, Poly$border[ind]),
               img    = Poly$img)
      }
    }
  }

  alpha <- dot$alpha; if (is.null(alpha)) alpha <- NA
  Poly$alpha   <- rep(alpha, length.out = length(Poly$col))
  if (! dot$shade$type == "none") 
    Poly <- color3D(Poly, plist$scalefac, dot$shade, lighting)

  Poly$proj   <- project(colMeans(Poly$x, na.rm = TRUE), 
                         colMeans(Poly$y, na.rm = TRUE), 
                         colMeans(Poly$z, na.rm = TRUE), plist)
      
  Poly$lwd     <- rep(lwd , length.out = length(Poly$col))
  Poly$lty     <- rep(lty , length.out = length(Poly$col))
  Poly$isimg   <- rep(0 , length.out = length(Poly$col))
  class(Poly)  <- "poly"

  if (image$add) 
    Poly <- XYimage (Poly, image, x, y, z, plist, col, breaks = breaks)

  if (contour$add) 
    segm <- contourfunc(contour, x, y, z, plist, cv, clim)
  else
    segm <- NULL

  if (iscolkey) 
    plist <- plistcolkey(plist, colkey, col, clim, clab, dot$clog, 
      type = "ribbon3D", breaks = breaks)
   
  plist <- plot.struct.3D(plist, poly = Poly, segm = segm, plot = plot)  

  setplist(plist)
  invisible(plist$mat)
}
