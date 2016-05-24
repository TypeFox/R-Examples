lines2D <- function(x, y, ...) {
  dot <- list(...)
  if (is.null(dot$type)) 
    dot$type <- "l"
  do.call("scatter2D", c(alist(x, y), dot))
}

points2D <- function(x, y,  ...) {
  dot <- list(...)
  if (is.null(dot$type)) 
    dot$type <- "p"
  do.call("scatter2D", c(alist(x, y), dot))
}


## =============================================================================
## Scatterplots (2-D)
## =============================================================================

scatter2D <- function(x, y, ..., colvar = NULL, 
                    col = NULL, NAcol = "white", breaks = NULL,
                    colkey = NULL, 
                    clim = NULL, clab = NULL, CI = NULL, 
                    add = FALSE, plot = TRUE) {

  plist <- initplist(add)

  plist <- add2Dplist(plist, "scatter", x = x, y = y, colvar = colvar, 
                    col = col, NAcol = NAcol, breaks = breaks,
                    colkey = colkey, 
                    clim = clim, clab = clab, CI = CI, ...)
  setplist(plist)
  if (!plot) return()
    
  dots <- splitpardots(list(...))

  isCI <- is.list(CI)
  if (isCI) 
    CI <- check.CI(CI, length(x), 2)


  if (is.null(col) & is.null(breaks))
     col <- jet.col(100)
   else if (is.null(col))
     col <- jet.col(length(breaks)-1)
  breaks <- check.breaks(breaks, col)
  
  if (! is.null(colvar)) {

    if (dots$clog) {
      colvar <- log(colvar)
      if (! is.null(clim)) clim <- log(clim) 
    }
    
    iscolkey <- is.colkey(colkey, col)        

    if (iscolkey) {
      colkey <- check.colkey(colkey)
      if (! add)       
        plist$plt$main <- colkey$parplt
      setplist(plist)    
      colkey$breaks <- breaks
    }  

    if (length(colvar) != length(x)) 
      stop ("length of 'colvar' should be equal to length of 'x' and 'y'")

    if (is.null(clim)) 
      clim <- range(colvar, na.rm = TRUE)
    
    if (! is.null(dots$alpha)) 
      col <- setalpha(col, dots$alpha)
    
    Col <- variablecol(colvar, col, NAcol, clim, breaks)

  } else  {  # no colvar
    Col <- col
    if (is.null(Col)) 
      Col <- "black"
    if (! is.null(dots$alpha)) 
      Col <- setalpha(Col, dots$alpha)
    iscolkey <- FALSE
  }   

  useSegments <- FALSE
  par (plt = plist$plt$main)

  if (! is.null(dots$points$type))
    if (dots$points$type %in% c("b", "l", "o"))
      if (length(Col) > 1  )
        useSegments <- TRUE

  if (useSegments) {
    Type <- dots$points$type
    len <- length(x)
    if (Type %in% c("b", "o"))   # no distinction is made..
      dots$points$type <- "p" 
    else     
      dots$points$type <- "n" 

   # mean of point colors for line colors
    LCol <- cbind(Col[-1], Col[-len])
    LCol <- apply(LCol, MARGIN = 1, FUN = MeanColors)
    if (! is.null(dots$alpha)) 
      LCol <- setalpha(LCol, dots$alpha)

    if (! add) 
      dots$main <- start2Dplot(dots$main, x, y)
    add <- TRUE

    if (isCI) {
      plot.CI.2d(CI, x, y, Col) 
      isCI <- FALSE
    }
    
    do.call("points", c(alist(x, y, col = Col), dots$points)) 
    dots$points$type <- NULL
    
    do.call("segments", c(alist(x[-len], y[-len], x[-1], y[-1], 
                  col = LCol), dots$points))
  }
  
  else if (! add) {
    dots$main <- start2Dplot(dots$main, x, y)
    
    if (isCI) {
      plot.CI.2d(CI, x, y, Col)   
      isCI <- FALSE
    }
    do.call("points", c(alist(x, y, col = Col), dots$points))
  } else  {
    
    if (isCI) {
      plot.CI.2d(CI, x, y, Col) 
      isCI <- FALSE
    }
    do.call("points", c(alist(x, y, col = Col), dots$points))
  }
    
  if (iscolkey) {
    drawcolkey(colkey, col, clim, clab, dots$clog) 
    par(plt = plist$plt$ori)  
  }    
  par(mar = par("mar"))

}

## =============================================================================
## Confidence interval check for scatters (2D and 3D)
## =============================================================================

check.CI <- function(CI, len, dim) {
  
    if (dim == 2 & !is.null(CI$z))
      stop("'CI' should not contain confidence intervals in the 'z' direction for 2-D plot")

    if (dim == 2 & is.null(CI$x) & is.null(CI$y))
      stop("'CI' should contain confidence intervals the 'x' and/or 'y' direction")
    
    if (dim == 3 & (is.null(CI$x) & is.null(CI$y) & is.null(CI$z)))
        stop("'CI' should contain confidence intervals in 'x', 'y', or 'z' direction")

    if (!is.null(CI$x)) {
      if (! is.matrix(CI$x))
        stop("'CI$x' should be a matrix")
      if (ncol(CI$x) != 2)
        stop("'CI$x' should be a matrix with two columns, with lower and upper value")
      if (nrow(CI$x) != len)
        stop("number of rows of matrix 'CI$x' should be equal to number of points")
    }      
    if (!is.null(CI$y)) {
      if (! is.matrix(CI$y))
        stop("'CI$y' should be a matrix")
      if (ncol(CI$y) != 2)
        stop("'CI$y' should be a matrix with two columns, with lower and upper value")
      if (nrow(CI$y) != len)
        stop("number of rows of matrix 'CI$y' should be equal to number of points")
    }      
    if (!is.null(CI$z)) {
      if (! is.matrix(CI$z))
        stop("'CI$z' should be a matrix")
      if (ncol(CI$z) != 2)
        stop("'CI$z' should be a matrix with two columns, with lower and upper value")
      if (nrow(CI$z) != len)
        stop("number of rows of matrix 'CI$z' should be equal to number of points")
    }      
  parameter <- list(alen = 0.01, lty = par("lty"), lwd = par("lwd"), col = NULL)
  CIpar <- CI
  CIpar$x <- CIpar$y <- CIpar$z <- NULL                           
  
  CIpar <- overrulepar(parameter, CIpar)
  CIpar$x <- CI$x
  CIpar$y <- CI$y
  CIpar$z <- CI$z
  CIpar
}


## =============================================================================
## CI in 2-d
## =============================================================================

plot.CI.2d <- function(CI, x, y, Col) {         # very-very simple

  CIpar <- CI[c("lty", "lwd", "col")]

  if (is.null(CIpar$col))
    CIpar$col <- Col
  
  if (! is.null(CI$x)) { # CI in x direction
     len <-  par("fin")[1] * CI$alen
     do.call("arrows", c(alist(x, y, x-CI$x[, 1], y, angle = 90, length = len), CIpar))
     do.call("arrows", c(alist(x, y, x+CI$x[, 2], y, angle = 90, length = len), CIpar))
  }  
  if (! is.null(CI$y)) { 
     len <-  par("fin")[2] * CI$alen
     do.call("arrows", c(alist(x, y, x, y-CI$y[, 1], angle = 90, length = len), CIpar))
     do.call("arrows", c(alist(x, y, x, y+CI$y[, 2], angle = 90, length = len), CIpar))
  }  
}

