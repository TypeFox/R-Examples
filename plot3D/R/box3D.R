
## =============================================================================
## box3D, rect3D, cube3D
## =============================================================================

## =============================================================================
## 3-D box function
## =============================================================================

border3D  <- function(x0, y0, z0, x1, y1, z1,
                   ..., colvar = NULL, phi = 40, theta = 40,
                   col = NULL, NAcol = "white", breaks = NULL,
                   colkey = NULL, 
                   panel.first = NULL,
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

  if (is.null(col) & is.null(breaks))
    col <- jet.col(100)
  else if (is.null(col))
    col <- jet.col(length(breaks)-1)

  breaks <- check.breaks(breaks, col)


  if (ispresent(colvar)) {
    if (length(colvar) != len)
      stop("'colvar' should have same length as 'x0', 'y0' and 'z0'")

    if (length(col) == 1)
      col <- c(col, col)

    if (is.null(clim))
      clim <- range(colvar, na.rm = TRUE)

    if (dot$clog) {      
      colvar <- log(colvar)
      clim <- log(clim)
    }

    iscolkey <- is.colkey(colkey, col)
    if (iscolkey) {
      colkey <- check.colkey(colkey)
      colkey$breaks <- breaks
    }

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
               phi = phi, theta = theta, plot = plot, col = col), dot$persp))
    plist <- getplist()
  }  

  if (is.function(panel.first))
    panel.first(plist$mat)

  lwd <- dot$points$lwd
  if (is.null(lwd)) 
    lwd <- 1
  lwd <- rep(lwd, length.out = len)

  lty <- dot$points$lty 
  if (is.null(lty)) 
    lty <- 1
  lty <- rep(lty, length.out = len)
  alpha <- dot$alpha; if (is.null(alpha)) alpha <- NA
  alpha <- rep(alpha, length.out = len)

  segm <- list()
  
 # small enough to use a loop 
  for (i in 1: len) {
    segm <- list(
      x.from = c(segm$x.from, x0[i], x1[i], x1[i], x0[i], x0[i],
        x1[i], x1[i], x0[i], x0[i], x1[i], x1[i], x0[i]),
      x.to   = c(segm$x.to  , x1[i], x1[i], x0[i], x0[i], x0[i],
        x1[i], x1[i], x0[i], x1[i], x1[i], x0[i], x0[i]),
      y.from = c(segm$y.from, y0[i], y0[i], y1[i], y1[i], y0[i],
        y0[i], y1[i], y1[i], y0[i], y0[i], y1[i], y1[i]),
      y.to   = c(segm$y.to  , y0[i], y1[i], y1[i], y0[i], y0[i],
        y0[i], y1[i], y1[i], y0[i], y1[i], y1[i], y0[i]),
      z.from = c(segm$z.from, z0[i], z0[i], z0[i], z0[i], z0[i],
        z0[i], z0[i], z0[i], z1[i], z1[i], z1[i], z1[i]),
      z.to   = c(segm$z.to  , z0[i], z0[i], z0[i], z0[i], z1[i],
        z1[i], z1[i], z1[i], z1[i], z1[i], z1[i], z1[i]),
      col    = c(segm$col, rep(Col[i], 12)),
      lwd    = c(segm$lwd, rep(lwd[i], length.out = 12)),
      lty    = c(segm$lty, rep(lty[i], length.out = 12)),
      alpha  = c(segm$alpha, alpha)
    )
  }
  
 # sort according to view
  segm$proj <- project(0.5*(segm$x.from + segm$x.to),
                       0.5*(segm$y.from + segm$y.to), 
                       0.5*(segm$z.from + segm$z.to), plist)
                  
  class(segm) <- "segments"

  if (iscolkey) 
    plist <- plistcolkey(plist, colkey, col, clim, clab, 
      dot$clog, type = "border3D", breaks = breaks)

  plist <- plot.struct.3D(plist, segm = segm, plot = plot)

  setplist(plist)   
  invisible(plist$mat)
}

## =============================================================================
## 3-D rect function
## =============================================================================

rect3D  <- function(x0, y0, z0, x1 = NULL, y1 = NULL, z1 = NULL,
                   ..., colvar = NULL, phi = 40, theta = 40,
                   col = NULL, NAcol = "white", breaks = NULL,
                   border = NA, facets = TRUE,
                   colkey = NULL, panel.first = NULL,
                   clim = NULL, clab = NULL, bty = "b", 
                   add = FALSE, plot = TRUE)  {

  plist <- initplist(add)

  dot  <- splitdotpersp(list(...), bty, NULL, c(x0, x1),
    c(y0, y1), c(z0, z1), plist = plist, breaks = breaks)

  len <- length(x0)
  if (length(y0) != len)
    stop("'y0' should have same length as 'x0'")
  if (length(z0) != len)
    stop("'z0' should have same length as 'x0'")

  numNULL <- is.null(x1) +  is.null(y1) + is.null(z1)
  if (numNULL != 1)
    stop ("one of 'x1', 'y1', or 'z1' should be NULL, and there are ", numNULL)
  if (! is.null(x1))
    if (length(x1) != len)
      stop("'x1' should have same length as 'x0' if not NULL")

  if (! is.null(y1))
    if (length(y1) != len)
      stop("'y1' should have same length as 'x0' if not NULL")

  if (! is.null(z1))
    if (length(z1) != len)
      stop("'z1' should have same length as 'x0' if not NULL")

  if (is.null(col) & is.null(breaks))
    col <- jet.col(100)
  else if (is.null(col))
    col <- jet.col(length(breaks)-1)

  breaks <- check.breaks(breaks, col)

  if (ispresent(colvar)) {
    if (length(colvar) != len)
      stop("'colvar' should have same length as 'x0', 'y0' and 'z0'")

    if (length(col) == 1)
      col <- c(col, col)

    if (is.null(clim))
      clim <- range(colvar, na.rm = TRUE)

    if (dot$clog) {                    # log transformation of color-values
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
      col <- "grey"
    if (! is.null(dot$alpha)) 
      col <- setalpha(col, dot$alpha)
    Col <- rep(col, length.out = len)
    iscolkey <- FALSE
  }

  Col <- createcolors(facets, border, Col)  

  if (is.null(plist)) {
    do.call("perspbox",
       c(alist(x = range(c(x0, x1)), y = range(c(y0, y1)),
               z = range(c(z0, z1)),
               phi = phi, theta = theta, plot = plot, col = col), dot$persp))
    plist <- getplist()
  }
    
  if (is.function(panel.first))
    panel.first(plist$mat)

  lwd <- dot$points$lwd
  if (is.null(lwd)) 
    lwd <- 1
  lwd <- rep(lwd, length.out = len)
  
  lty <- dot$points$lty 
  if (is.null(lty)) 
    lty <- 1
  lty <- rep(lty, length.out = len)

  alpha <- dot$alpha; if (is.null(alpha)) alpha <- NA
  alpha <- rep(alpha, length.out = len)

  poly <- list()
  for (i in 1: len) {
    if (is.null(x1)) {
      x = c(x0[i], x0[i], x0[i], x0[i], NA)
      y = c(y0[i], y1[i], y1[i], y0[i], NA)
      z = c(z0[i], z0[i], z1[i], z1[i], NA)
    } else if (is.null(y1)) {  
      x = c(x0[i], x1[i], x1[i], x0[i], NA)
      y = c(y0[i], y0[i], y0[i], y0[i], NA)
      z = c(z0[i], z0[i], z1[i], z1[i], NA)
    } else if (is.null(z1)) {  
      x = c(x0[i], x0[i], x1[i], x1[i], NA)
      y = c(y0[i], y1[i], y1[i], y0[i], NA)
      z = c(z0[i], z0[i], z0[i], z0[i], NA)
    } 
    poly <- list(
      x      = cbind(poly$x, x),
      y      = cbind(poly$y, y),
      z      = cbind(poly$z, z),
      col    = c(poly$col,    Col$facet[i]),
      border = c(poly$border, Col$border[i]),
      lwd    = c(poly$lwd,    lwd[i]),
      lty    = c(poly$lty,    lty[i]),
      isimg  = c(poly$isimg, 0),
      alpha  = c(poly$alpha, alpha)
      
    )
  }
 # projection depth
  poly$proj <- project(colMeans(poly$x, na.rm = TRUE),
                       colMeans(poly$y, na.rm = TRUE), 
                       colMeans(poly$z, na.rm = TRUE), plist, FALSE)
  
  class(poly) <- "poly"

  if (iscolkey ) 
    plist <- plistcolkey(plist, colkey, col, clim, clab, 
      dot$clog, type = "rect3D", breaks = breaks)

  plist <- plot.struct.3D(plist, poly = poly, plot = plot)

  setplist(plist)
  invisible(plist$mat)
}

## =============================================================================
## 3-D box function with colored facets
## =============================================================================

box3D  <- function(x0, y0, z0, x1, y1, z1,
                   ...,  colvar = NULL, phi = 40, theta = 40,
                   col = NULL, NAcol = "white", breaks = NULL,
                   border = NA, facets = TRUE,
                   colkey = NULL, panel.first = NULL,
                   clim = NULL, clab = NULL, bty = "b", 
                   add = FALSE, plot = TRUE)  {

  plist <- initplist(add)

  dot  <- splitdotpersp(list(...), bty, NULL, c(x0, x1),
    c(y0, y1), c(z0, z1), plist = plist, breaks = breaks)

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

  if (is.null(col) & is.null(breaks))
    col <- jet.col(100)
  else if (is.null(col))
    col <- jet.col(length(breaks)-1)

  breaks <- check.breaks(breaks, col)

  if (ispresent(colvar)) {
    if (length(colvar) != len)
      stop("'colvar' should have same length as 'x0', 'y0' and 'z0'")

    if (length(col) == 1)
      col <- c(col, col)

    if (is.null(clim))
      clim <- range(colvar, na.rm = TRUE)

    if (dot$clog) {                    # log transformation of color-values
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
      col <- "grey"
    if (! is.null(dot$alpha)) 
      col <- setalpha(col, dot$alpha)
    Col <- rep(col, length.out = len)
    iscolkey <- FALSE
  }

# The colors of facets and border
  Col <- createcolors(facets, border, Col)  

  if (is.null(plist)) {
    do.call("perspbox",
       c(alist(x = range(c(x0, x1)), y = range(c(y0, y1)),
               z = range(c(z0, z1)),
               phi = phi, theta = theta, plot = plot, col = col), dot$persp))
    plist <- getplist()
  }
    
  if (is.function(panel.first))
    panel.first(plist$mat)

  lwd <- dot$points$lwd
  if (is.null(lwd)) 
    lwd <- 1
  lwd <- rep(lwd, length.out = len)

  lty <- dot$points$lty 
  if (is.null(lty)) 
    lty <- 1
  lty <- rep(lty, length.out = len)

  alpha <- dot$alpha; if (is.null(alpha)) alpha <- NA
  alpha <- rep(alpha, length.out = len)

  poly <- list()
  for (i in 1: len) {
    poly <- list(
      x = cbind(poly$x, c(x0[i], x1[i], x1[i], x0[i], NA), 
                        c(x1[i], x1[i], x1[i], x1[i], NA), 
                        c(x1[i], x0[i], x0[i], x1[i], NA),
                        c(x0[i], x0[i], x0[i], x0[i], NA),
                        c(x0[i], x1[i], x1[i], x0[i], NA),
                        c(x0[i], x1[i], x1[i], x0[i], NA)),
      y = cbind(poly$y, c(y0[i], y0[i], y0[i], y0[i], NA), 
                        c(y0[i], y1[i], y1[i], y0[i], NA),
                        c(y1[i], y1[i], y1[i], y1[i], NA), 
                        c(y1[i], y0[i], y0[i], y1[i], NA),
                        c(y0[i], y0[i], y1[i], y1[i], NA),
                        c(y0[i], y0[i], y1[i], y1[i], NA)),
      z = cbind(poly$z, c(z0[i], z0[i], z1[i], z1[i], NA),
                        c(z0[i], z0[i], z1[i], z1[i], NA),
                        c(z0[i], z0[i], z1[i], z1[i], NA), 
                        c(z0[i], z0[i], z1[i], z1[i], NA),
                        c(z1[i], z1[i], z1[i], z1[i], NA), 
                        c(z0[i], z0[i], z0[i], z0[i], NA)),
      col    = c(poly$col,    rep(Col$facet[i], length.out = 6)),
      border = c(poly$border, rep(Col$border[i], length.out = 6)),
      lwd    = c(poly$lwd,    rep(lwd[i], length.out = 6)),
      lty    = c(poly$lty,    rep(lty[i], length.out = 6)),
      isimg  = c(poly$isimg,  rep(0, length.out = 6)),
      alpha  = c(poly$alpha, rep(alpha, length.out = 6))
      
    )
  }
 # sort according to view
  poly$proj <- project(colMeans(poly$x, na.rm = TRUE),
                       colMeans(poly$y, na.rm = TRUE), 
                       colMeans(poly$z, na.rm = TRUE), plist, FALSE)
  class(poly) <- "poly"

  if (iscolkey ) 
    plist <- plistcolkey(plist, colkey, col, clim, clab, 
      dot$clog, type = "box3D", breaks = breaks)

  plist <- plot.struct.3D(plist, poly = poly, plot = plot)

  setplist(plist)
  invisible(plist$mat)
}











