## =============================================================================
## 3-d spherical surface
## =============================================================================

spheresurf3D <- function(colvar = matrix(nrow = 50, ncol = 50, data = 1:50, byrow = TRUE),
     ..., phi = 0, theta = 0,
     col = NULL, NAcol = "white", breaks = NULL,
     border = NA, facets = TRUE,
     contour = FALSE, colkey = NULL, resfac = 1,
     panel.first = NULL, clim = NULL, clab = NULL,
     bty = "n", lighting = FALSE, shade = NA, ltheta = -135, lphi = 0,
     inttype = 1, full = FALSE, add = FALSE, plot = TRUE) {

  plist <- initplist(add)

  r <- 1  # sphere radius
  
  if (!is.matrix(colvar))
    stop("'colvar' should be a matrix or absent")

  X <- seq(0, 2*pi, length = nrow(colvar))
  Y <- seq(0,   pi, length = ncol(colvar))
  
 # change resolution
  if (any(resfac != 1)) {
    res <- changeres(resfac, X, Y, colvar)
    X <- res$x
    Y <- res$y
    colvar <- res$z
  }

 # check contours
  contour <- check.args(contour)
  
  cv <- colvar

  M <- mesh(X, Y)
  x <- with (M, -r*cos(x)*sin(y))
  y <- with (M, -r*sin(y)*sin(x))
  z <- with (M, -r*cos(y))
  
  dot <- splitdotpersp(list(...), bty, lighting, 
    x, y, z, plist = plist, shade, lphi, ltheta, breaks = breaks)

  DD <- dim(x)

  if (is.null(col))
    if (is.null(breaks))
      col <- jet.col(100)
    else
      col <- jet.col(length(breaks)-1)

  CC <- check.colvar.persp(colvar, z, col, inttype, clim, dot$alpha)
  colvar <- CC$colvar
  col <- CC$col
  
  if (is.null(clim))
    clim <- range(colvar, na.rm = TRUE)
  
  if (dot$clog) {
    colvar <- log(colvar)
    clim <- log(clim)
  }

  iscolkey <- is.colkey(colkey, col)
  if (iscolkey) 
    colkey <- check.colkey(colkey)

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

  Poly <- paintit(colvar, x, y, z, plist, col, NAcol, clim, border,
          facets, dot$points$lwd, dot$points$lty, dot, Extend, !full, 
          breaks = breaks)

  if (contour$add) {
    contour$side <- NULL #"z"
    Pmin <- min(Poly$img[[1]]$sl$Proj) # minimal projection depth behind which lines not drawn
    
    col.lines <- contour$args$col
    if (is.null(col.lines))
      col.lines <- "black"
    
    segm <- NULL
    if (is.null(    contour$args$nlevels))
      contour$args$nlevels <- 10
    if (is.null(    contour$args$levels))
      contour$args$levels <- pretty(range(cv, na.rm = TRUE), contour$args$nlevels)
    line.list <- do.call("contourLines",
      alist(x = X, y = Y, z = cv, nlevels = contour$args$nlevels, 
      levels = contour$args$levels))

    contour$args$nlevels <- contour$args$col <- contour$args$levels <- NULL

    for (i in 1:length(line.list)) {
       clines <- line.list[[i]]
       X <- clines$x; Y <- clines$y
       x <-  -r*cos(X)*sin(Y)
       y <-  -r*sin(Y)*sin(X)
       z <-  -r*cos(Y)
       sl <- sortlistvec(x, y, z, plist, ignorez = FALSE)

       isel <- which(sl$Proj > Pmin)

       if (length(isel) > 1) {
        segm <- do.call("addlines", c(alist(segm, x[isel], y[isel], 
          z[isel], col = col.lines, plist = plist, ignorez = FALSE), contour$args))
       }
    }
    if (! is.null(segm))
      segm$proj <- segm$proj + 1e-1  # put it on foreground...

  } else
    segm <- NULL

  if (iscolkey) 
    plist <- plistcolkey(plist, colkey, col, clim, clab, 
         dot$clog, type = "spheresurf3D", breaks = breaks)

  plist <- plot.struct.3D(plist, poly = Poly, segm = segm, plot = plot)  

  setplist(plist)   
  invisible(plist$mat)
}
