                                       
# ==============================================================================
# the 2D image and contour function, using rgl
# ==============================================================================

image2Drgl <- function(z, x = seq(0, 1, length.out = nrow(z)),
              y = seq(0, 1, length.out = ncol(z)), ...,
              col = NULL, NAcol = "white", breaks = NULL,
              border = NA, facets = TRUE, colkey = NULL,
              contour = FALSE, smooth = FALSE,
              clim = NULL, clab = NULL, shade = NA,
              inttype = 1, dz = 0, add = FALSE)  {
# ------------------------------------------------------------------------------

  if (! is.null(list(...)$rasterImage))
    smooth <- list(...)$rasterImage

  namesextra <- c("border", "inttype", "facets")
  if (is.null(col) & is.null(breaks))
    col <- jet.col(100)
  else if (is.null(col))
    col <- jet.col(length(breaks)-1)
  zz <- 1 + dz
  if (!is.na(shade)) {
  dots <- plot2Drgl("persp3D", x, y, z, colvar = z, col = col, NAcol = NAcol,
     breaks = breaks, clim = clim, add = add, clab = clab,
     namesextra = namesextra,  shade = shade, border = border,
     colkey = colkey, inttype = inttype, smooth = smooth, ...)
  zz <- z
  }
  else
  dots <- plot2Drgl("image3D", x, y, colvar = z, col, NAcol,
     breaks = breaks, clim, add, clab = clab,
     namesextra = namesextra,  z = zz, border = border,
     colkey = colkey, inttype = inttype, smooth = smooth, ...)
                    
  iscontour <- contour
  if (! is.logical(iscontour))
    iscontour <- TRUE
  else
    contour <- list()
  if (iscontour)  {
    if (is.null(contour$col))
      contour$col <- "black"
    do.call("contour3D", c(alist(z = zz + 0.001, x = x, y = y, colvar = z,
      add = TRUE, plot = FALSE), contour, colkey = FALSE))
  }

  finishplotrgl(dots, namesextra = namesextra)
  
}

# ==============================================================================

contour2Drgl <- function(z, x = seq(0, 1, length.out = nrow(z)),
                    y = seq(0, 1,  length.out = ncol(z)), ...,
                    col = NULL, colkey = NULL,  clim = NULL,
                    clab = NULL, dz = 0.1, add = FALSE)  {
# ------------------------------------------------------------------------------

  namesextra <- c("lty", "lwd", "nlevels", "levels")
  if (! is.null(list(...)$breaks))
    stop("'breaks' cannot be used with contour2Drgl - use levels instead")
  dots <- plot2Drgl("contour3D", x, y, z = 1 + dz,
     colvar = z, col = col, breaks = NULL, colkey = colkey, clim = clim,
     clab = clab, namesextra = namesextra, add = add, ...)
  finishplotrgl(dots, namesextra = namesextra)
}



