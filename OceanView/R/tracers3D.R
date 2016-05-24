
## =============================================================================
## =============================================================================
## 3D tracer distributions
## =============================================================================
## =============================================================================

## =============================================================================
## Traditional graphics
## =============================================================================

tracers3D <- function(x, y, z, colvar = NULL, ..., 
    col = NULL, NAcol = "white", breaks = NULL,
    colkey = FALSE, clim = NULL, clab = NULL, surf = NULL) {

  if (!is.null(surf))  {
    surf$plot <- FALSE
    surf$add <- FALSE
    do.call("persp3D", surf)
  }
  
  plist <- getplist()

  dots <- list(...)

  if (!is.null(dots$main))  {
    plist$dot$main <- dots$main
    dots$main <- NULL
  }
  
  plist$pt <- NULL
 
  if (! is.null(plist$numkeys))
    if (plist$numkeys > 0) {
      for (i in plist$numkeys:1)
        if (plist$colkey[[i]]$type == "scatter3D") 
           plist$colkey[[i]] <- NULL
      plist$numkeys <- length(plist$colkey)
    }
  
  setplist(plist)
  do.call("points3D", c(alist(x, y, z, colvar = colvar, 
    col = col, NAcol = NAcol, breaks = breaks, clim = clim, clab = clab,
    add = TRUE, plot = FALSE, colkey = colkey), dots))
  plotdev()
}

## =============================================================================
## same in Open GL graphics
## =============================================================================

tracers3Drgl <- function(x, y, z, colvar = NULL, ..., 
  col = NULL, NAcol = "white", breaks = NULL,
  colkey = FALSE, clim = NULL, clab = NULL) {

  x <- as.vector(x)
  y <- as.vector(y)
  z <- as.vector(z)
  
  len <- length(x)
  if (length(y) != len)
    stop("'y' should be of same length as 'x'")
  if (length(z) != len)
    stop("'z' should be of same length as 'x'")
    
  save <- par3d(skipRedraw = TRUE, ignoreExtent = TRUE)
  on.exit(par3d(save))
  
  dots <- list(...)

  if (is.null(col))
    col <- "black"
  breaks <- check.breaks(breaks, col)

  # colors and color variable
  if (! is.null(colvar)) {
    if (length(colvar) != len) 
      stop ("dimension of 'colvar' should be equal to dimension of 'x', 'y' and 'z'")

    if (is.null(clim)) 
      clim <- range(colvar, na.rm = TRUE)
    
    Col <- variablecol(colvar, col, NAcol, clim, breaks) # generate color scheme

  } else {
    Col <- col
  }   

  cex <- dots$cex
  if (is.null(cex )) cex <- 1
  
  alpha <- dots$alpha
  if (is.null(alpha)) 
    alpha <- material3d()$alpha

 # if main is passed...
  plist <- getplist()
  if (!is.null(dots$main))  {
    if (is.null(plist$dot$main))
      plist$dot$main <- dots$main
    ids <- plist$rgl$D
    if (! is.null(ids))
      if (!is.null(ids$main)) 
        rgl.pop(type = "shapes", id = ids["main"])
    M <- mtext3d(dots$main, "x++", line = 2)
    dots$main <- NULL
    plist$rgl$D$main <- M
  }

  pp <- rgl.ids()
#  rgl.pop(type = "shapes", id = pp[which(as.character(pp$type) == "spheres"),1])
#  spheres3d(x, y, z, radius = cex * 0.0175 * plist$scalefac$expand, col = Col)
  rgl.pop(type = "shapes", id = pp[which(as.character(pp$type) == "points"),1])
  plot3d(x = x, y = y, z = z,
               size = 6 *cex, col = Col, add = TRUE, alpha = alpha)
               
  plist$pt <- list(x.mid = x, y.mid = y, z.mid = z,
    col = Col, pch = rep(1, length(x)), bg = rep(1, length(x)), 
    cex = rep(cex, length(x)), alpha = rep(alpha, length(x)),
    proj = rep(NA, length(x)))
   
  setplist(plist)
}


