
## =============================================================================
## text in 3D
## =============================================================================
# x, y, z, colvar: same length
    
text3D <- function(x, y, z, labels, ..., colvar = NULL, 
                   phi = 40, theta = 40,
                   col = NULL, NAcol = "white", breaks = NULL,
                   colkey = NULL, 
                   panel.first = NULL,
                   clim = NULL, clab = NULL, bty = "b",
                   add = FALSE, plot = TRUE) {

  plist <- initplist(add)

  x <- as.vector(x)
  y <- as.vector(y)
  z <- as.vector(z)
  
  if (length(y) != length(x))
    stop("'y' should have same length as 'x'")
  if (length(z) != length(x))
    stop("'z' should have same length as 'x'")
  if (length(labels) != length(x))
    stop("'labels' should have same length as 'x'")

  dot <- splitdotpersp(list(...), bty, NULL, x, y, z, plist = plist, breaks = breaks)

  if (ispresent(colvar)) {
  
    if (length(colvar) != length(x))
      stop("'colvar' should have same length as 'x', 'y' and 'z'")

    colvar <- as.vector(colvar)
    
    if (is.null(clim)) 
      clim <- range(colvar, na.rm = TRUE)
    
    if (dot$clog) {
      colvar <- log(colvar)
      clim <- log(clim) 
    }

    if (is.null(col))
      if (is.null(breaks))
        col <- jet.col(100)
      else
        col <- jet.col(length(breaks)-1)

    breaks <- check.breaks(breaks, col)
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
    Col <- rep(col, length.out = length(x))  
    iscolkey <- FALSE
  }
  
  if (is.null(plist)) {
    do.call("perspbox", c(alist(x = range(x), y = range(y), 
             z = range(z, na.rm = TRUE),
             phi = phi, theta = theta, plot = plot, colkey = colkey, col = col), 
             dot$persp))
    plist <- getplist()
  }
    
  if (is.function(panel.first)) 
    panel.first(plist$mat)  

  Proj   <- project(x, y, z, plist)      # sort labels according to view

  setargs <- function(dot, default) {
    if (is.null(dot)) 
      rep(default, length.out = length(x))
    else if (is.vector(dot) & length(dot) > 1)
      stop("cannot use a vector for arguments of 'text3D'")
    else 
      rep(unlist(dot), length.out = length(x))  
  }
  alpha <- dot$alpha; if (is.null(alpha)) alpha <- NA
  alpha <- rep(alpha, length.out = length(x))

  tlist <- list(x    = x,
                y    = y,
                z    = z,                                  
                labels = labels,
                col  = Col,
                adj = setargs (dot$points$adj, 0),
                cex = setargs (dot$points$cex, 1),
                font = setargs(dot$points$font, 1),
                srt = setargs(dot$points$srt, 0),
                alpha = alpha,
                proj = Proj)                 

  if (iscolkey) 
    plist <- plistcolkey(plist, colkey, col, clim, clab, 
         dot$clog, type = "label3D", breaks = breaks)
                 
  plist <- plot.struct.3D (plist, labels = tlist, plot = plot)        

  setplist(plist)   
  invisible(plist$mat)
}





## =============================================================================
## text in 2D
## =============================================================================
    
text2D <- function(x, y, labels, ..., colvar = NULL, 
                   col = NULL, NAcol = "white", breaks = NULL,
                   colkey = NULL, 
                   clim = NULL, clab = NULL, add = FALSE, plot = TRUE) {

  if (add) 
    plist <- getplist()
  else
    plist <- NULL

  plist <- add2Dplist(plist, "text", x = x, y = y, labels = labels,
                    colvar = colvar, col = col, NAcol = NAcol,
                    breaks = breaks, colkey = colkey,
                    clim = clim, clab = clab, ...)
  setplist(plist)
  if (!plot) return()

  x     <- as.vector(x)
  y     <- as.vector(y)
  
  if (length(y) != length(x))
    stop("'y' should have same length as 'x'")

  if (length(labels) != length(x))
    stop("'labels' should have same length as 'x'")

  dots <- splitpardots(list(...))

  if (! is.null(colvar)) {
    if (is.null(col))
      if (is.null(breaks))
        col <- jet.col(100)
      else
        col <- jet.col(length(breaks)-1)

    if (dots$clog) {
      colvar <- log(colvar)
      if (! is.null(clim)) clim <- log(clim)
    }

    iscolkey <- is.colkey(colkey, col)

    if (iscolkey) {
      colkey <- check.colkey(colkey, add)
      if (! add)
        par.ori <- par(plt = colkey$parplt)
      colkey$breaks <- breaks
    }

    if (length(colvar) != length(x))
      stop ("length of 'colvar' should be equal to length of 'x', and 'y'")

    if (is.null(clim))
      clim <- range(colvar, na.rm = TRUE)

    if (! is.null(dots$alpha)) col <- setalpha(col, dots$alpha)
    Col <- variablecol(colvar, col, NAcol, clim, breaks)

  } else  {  
    Col <- col
    if (is.null(Col)) 
      Col <- "black"
    if (! is.null(dots$alpha)) 
      Col <- setalpha(Col, dots$alpha)
    iscolkey <- FALSE
  }

  if (! add)
    dots$main <- start2Dplot(dots$main, x, y)

  do.call("text", c(alist(x, y, labels = labels, col = Col), dots$points))

  if (iscolkey) {
    drawcolkey(colkey, col, clim, clab, dots$clog)
    if (! add)
      par(plt = par.ori)
    par(mar = par("mar"))
  }
}
