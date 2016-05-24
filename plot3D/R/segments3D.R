## =============================================================================
## 3-D segments function; no shading/light
## =============================================================================

segments3D  <- function(x0, y0, z0, x1 = x0, y1 = y0, z1 = z0,
                    ..., colvar = NULL, phi = 40, theta = 40,
                    col = NULL, NAcol = "white", breaks = NULL,
                    colkey = NULL, panel.first = NULL,
                    clim = NULL, clab = NULL, bty = "b", 
                    add = FALSE, plot = TRUE)  {

  plist <- initplist(add)

  dot  <- splitdotpersp(list(...), bty, NULL, 
    c(x0, x1), c(y0, y1), c(z0, z1), plist = plist, breaks = breaks)

  len <- length(x0)
  if (length(y0) != len)
    stop("'y0' should have same length as 'x0'")
  if (length(z0) != len)
    stop("'z0' should have same length as 'x0'")
  if (length(x1) != len)
    stop("'x1' should have same length as 'x0'")
  if (length(y1) != len)
    stop("'y1' should have same length as 'x0'")
  if (length(z1) != len)
    stop("'z1' should have same length as 'x0'")
  if (ispresent(colvar)) {
    if (length(colvar) != len)
      stop("'colvar' should have same length as 'x0', 'y0' and 'z0'")
    
    if (is.null(col) & is.null(breaks))
      col <- jet.col(100)
    else if (is.null(col))
      col <- jet.col(length(breaks)-1)
    breaks <- check.breaks(breaks, col)

    if (length(col) == 1)
      col <- c(col, col)

    if (is.null(clim)) 
      clim <- range(colvar, na.rm = TRUE)
    
    if (dot$clog) {                    
      colvar <- log(colvar)
      clim <- log(clim)
    }

    iscolkey <- is.colkey(colkey, col) 
    if (iscolkey) 
      colkey <- check.colkey(colkey)

    if (! is.null(dot$alpha)) 
      col <- setalpha(col, dot$alpha)
    Col <- variablecol(colvar, col, NAcol, clim, breaks)

  } else {
    if (is.null(col))
      col <- "black"
    if (! is.null(dot$alpha)) 
      col <- setalpha(col, dot$alpha)
    Col <- rep(col, length.out = len)  
    iscolkey <- FALSE
  }   

  if (is.null(plist)) {
    do.call("perspbox", 
       c(alist(x = range(c(x0, x1)), y = range(c(y0, y1)), 
               z = range(c(z0, z1)), 
               phi = phi, theta = theta, plot = plot, 
               colkey = colkey, col = col), dot$persp))
    plist <- getplist()
  }  
  if (is.function(panel.first)) 
    panel.first(plist$mat)
  
  lwd <- dot$points$lwd
  if (is.null(lwd)) 
    lwd <- 1

  lty <- dot$points$lty
  if (is.null(lty)) 
    lty <- 1

  alpha <- dot$alpha; if (is.null(alpha)) alpha <- NA
  alpha <- rep(alpha, length.out = len)

  Proj   <- project(0.5*(x0 + x1), 0.5*(y0 + y1), 
                    0.5*(z0 + z1), plist)

  segm <- list(x.from = x0, 
               x.to   = x1,
               y.from = y0, 
               y.to   = y1,                                  
               z.from = z0, 
               z.to   = z1,                                  
               col    = Col,
               lwd    = rep(lwd , length.out = len),
               lty    = rep(lty , length.out = len),
               alpha  = alpha,
               proj   = Proj)
  class(segm) <- "segments"

  if (iscolkey) 
    plist <- plistcolkey(plist, colkey, col, clim, clab, 
      dot$clog, type = "segments3D", breaks = breaks)

  plist <- plot.struct.3D(plist, segm = segm, plot = plot)  

  setplist(plist)                      
  invisible(plist$mat)
}

