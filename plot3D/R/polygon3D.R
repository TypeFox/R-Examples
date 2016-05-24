## =============================================================================
## 3-D polygon function
## =============================================================================

polygon3D  <- function(x, y, z, 
                    ..., colvar = NULL, phi = 40, theta = 40,
                    col = NULL, NAcol = "white", breaks = NULL,
                    border = NA, facets = TRUE,
                    colkey = NULL, panel.first = NULL,
                    clim = NULL, clab = NULL, bty = "b", 
                    add = FALSE, plot = TRUE)  {

  plist <- initplist(add)

  dot  <- splitdotpersp(list(...), bty, NULL, x, y, z, plist = plist, breaks = breaks)

  checkinput <- function (x) {
    if (is.matrix(x)) {
     x <- as.vector(x)
     if (is.na (x[length(x)]) )
       x <- x[-length(x)]
    }
    x
  }
  x <- checkinput(x)
  y <- checkinput(y)
  z <- checkinput(z)
  
  if (length(y) != length(x))
    stop("'y' should have same length as 'x'")
  if (length(z) != length(x))
    stop("'z' should have same length as 'x'")

 # check for NAs (and number of polygons)
  len <- 1
  if (any (is.na(x) | is.na(y) | is.na(z))) {
    i1 <- which(is.na(x))
    i2 <- which(is.na(y))
    i3 <- which(is.na(z))
    ii <- unique(c(i1, i2, i3))
    if (1 %in% ii | length(x) %in% ii)
      stop ("first or last element of 'x', 'y', or 'z' cannot be 'NA'")
    di <- diff(sort(c(0, ii, length(x)+1)))-1
    if (min(di) == 1)
      stop ("two consecutive elements of 'x', 'y', or 'z' cannot be 'NA'")
    x[ii] <- NA
    y[ii] <- NA
    z[ii] <- NA
    len <- length(ii) + 1  # number of polygons!

    xx <- yy <- zz <- matrix(nrow = max(di) + 1, ncol = len, data = NA)
    ii <- c(0, ii, length(x)+ 1)
    for (i in 1 : len) {
      iseq <- (ii[i]+1): (ii[i+1]-1)
      xx[1:length(iseq), i] <- x[iseq]
      yy[1:length(iseq), i] <- y[iseq]
      zz[1:length(iseq), i] <- z[iseq]
    }
  } else {
    xx <- matrix(ncol = 1, data = c(x, NA)) 
    yy <- matrix(ncol = 1, data = c(y, NA)) 
    zz <- matrix(ncol = 1, data = c(z, NA)) 
  }
  breaks <- check.breaks(breaks, col)

  if (ispresent(colvar)) { 
    if (length(colvar) != len)
      stop("'colvar' should have same length as number of polygons (= 1+ number of NAs in 'x', 'y' and 'z')")
    
    if (is.null(col))
      col <- jet.col(100)
    
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
           c(alist(x = range(x, na.rm = TRUE), y = range(y, na.rm = TRUE), 
             z = range(z, na.rm = TRUE), 
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

  Poly <- list(x = xx, 
               y = yy, 
               z = zz, 
               col    = Col$facet,
               border = Col$border,
               lwd    = rep(lwd , length.out = len),
               lty    = rep(lty , length.out = len),
               alpha  = alpha, 
               isimg  = rep(0, length.out = len))

  Poly$proj   <- project(colMeans(xx, na.rm = TRUE), colMeans(yy, na.rm = TRUE), 
    colMeans(zz, na.rm = TRUE), plist)

  class(Poly) <- "poly"

  if (iscolkey) 
    plist <- plistcolkey(plist, colkey, col, clim, clab, dot$clog, 
      type = "polygon3D", breaks = breaks)

  plist <- plot.struct.3D(plist, poly = Poly, plot = plot)  

  setplist(plist)   
  invisible(plist$mat)
}

## =============================================================================
## 2-D polygon and triangle function
## =============================================================================

polygon2D  <- function(x, y, ..., colvar = NULL, 
                    col = NULL, NAcol = "white", breaks = NULL,
                    border = NA, facets = TRUE,
                    colkey = NULL, 
                    clim = NULL, clab = NULL, add = FALSE, plot = TRUE)  {

  plist <- initplist(add)

  plist <- add2Dplist(plist, "polygon", x = x, y = y, colvar = colvar,
                    col = col, NAcol = NAcol, breaks = breaks,
                    border = border, facets = facets,
                    colkey = colkey, clim = clim,
                    clab = clab, ...)
  setplist(plist)
  if (! plot) return()
  
  dots <- splitpardots(list(...))

  checkinput <- function (x) {
    if (is.matrix(x)) {
     x <- as.vector(x)
     if (is.na (x[length(x)]) )
       x <- x[-length(x)]
    }
    x
  }
  x <- checkinput(x)
  y <- checkinput(y)
 
  if (length(y) != length(x))
    stop("'y' should have same length as 'x'")

 # check for NAs (and number of polygons)
  len <- 1
  if (any (is.na(x) | is.na(y))) {
    i1 <- which(is.na(x))
    i2 <- which(is.na(y))
    ii <- unique(c(i1, i2))
    if (1 %in% ii | length(x) %in% ii)
      stop ("first or last element of 'x', or 'y' cannot be 'NA'")
    di <- diff(sort(c(0, ii, length(x)+1)))-1
#    if (min(di) == 1)
#      stop ("two consecutive elements of 'x' or 'y' cannot be 'NA'")
    x[ii] <- NA
    y[ii] <- NA

    len <- length(ii) + 1  # number of polygons!
    
    xx <- yy <- matrix(nrow = max(di) + 1, ncol = len, data = NA)
    ii <- c(0, ii, length(x)+ 1)
    for (i in 1 : len) {
      iseq <- (ii[i]+1): (ii[i+1]-1)
      xx[1:length(iseq), i] <- x[iseq]
      yy[1:length(iseq), i] <- y[iseq]
    }
  } else {
    xx <- matrix(ncol = 1, data = c(x, NA)) 
    yy <- matrix(ncol = 1, data = c(y, NA)) 
  }
  breaks <- check.breaks(breaks, col)
  if (ispresent(colvar)) {
    if (length(colvar) != len)
      stop("'colvar' should have same length as number of polygons (= 1+ number of NAs in 'x', 'y' and 'z')")
    
    if (is.null(col))
      col <- jet.col(100)
    
    if (length(col) == 1)
      col <- c(col, col)

    if (is.null(clim)) 
      clim <- range(colvar, na.rm = TRUE)
    
    if (dots$clog) {                       
      colvar <- log(colvar)
      clim <- log(clim)
    }

    iscolkey <- is.colkey(colkey, col) 
    if (iscolkey) {
      colkey <- check.colkey(colkey, add)
      if (! add)
        par.ori <- par(plt = colkey$parplt)
      colkey$breaks <- breaks

    }
     
    if (! is.null(dots$alpha)) 
      col <- setalpha(col, dots$alpha)
    Col <- variablecol(colvar, col, NAcol, clim, breaks)

  } else {
    if (is.null(col))
      col <- "grey"
    if (! is.null(dots$alpha)) 
      col <- setalpha(col, dots$alpha)
    Col <- rep(col, length.out = len)  
    iscolkey <- FALSE
  }   

# The colors of facets and border
  Col <- createcolors(facets, border, Col)  

  if (! add)
    dots$main <- start2Dplot(dots$main, x, y)
    
  do.call("polygon", c(alist(x, y, col = Col$facet, border = Col$border), dots$points))

  if (iscolkey) {
    drawcolkey(colkey, col, clim, clab, dots$clog)
    if (! add)
      par(plt = par.ori)
    par(mar = par("mar"))
  }
}



