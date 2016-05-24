# ==============================================================================
# the 2D arrows function, using rgl
# ==============================================================================

arrows2Drgl <- function(x0, y0, x1, y1, colvar = NULL, ...,
              col = NULL, NAcol = "white", breaks = NULL, colkey = NULL,
              clim = NULL, clab = NULL, type = "simple", dz = 0.1, add = FALSE)  {
# ------------------------------------------------------------------------------

  namesextra <- c("code", "length", "angle", "lwd", "lty", "type")
  dots <- plot2Drglbis("arrows3D", x0 = x0, y0 = y0, x1 = x1, y1 = y1,
     dz = dz, colkey = colkey, type = type, colvar = colvar, col = col,
     NAcol = NAcol, breaks = breaks, clim = clim, clab = clab,
     add = add, namesextra = namesextra, ...)
  finishplotrgl(dots, namesextra)
}

# ==============================================================================

rect2Drgl <- function (x0, y0, x1, y1, colvar = NULL, ...,
         col = NULL, NAcol = "white", breaks = NULL,
         colkey = NULL, clim = NULL, clab = NULL,
         dz = 0.1, add = FALSE)   {


  namesextra <- c("lwd", "lty")
  z <- rep(1+dz, length.out = length(x0))
  dots <- plot2Drglbis("box3D", x0 = x0, y0 = y0, x1 = x1, y1 = y1,
    dz = dz, colkey = colkey, colvar = colvar, col = col,
    NAcol = NAcol, breaks = breaks, clim = clim, clab = clab, add = add,
    namesextra = namesextra, ...)
  finishplotrgl(dots, namesextra)
}

# ==============================================================================
# the 2D segments function, using rgl
# ==============================================================================

segments2Drgl <- function(x0, y0, x1, y1, colvar = NULL, ...,
               col = NULL, NAcol = "white", breaks = NULL, colkey = NULL,
               clim = NULL, clab = NULL, dz = 0.1, add = FALSE)  {
# ------------------------------------------------------------------------------

  namesextra <- c("lwd", "lty")
  z <- rep(1+dz, length.out = length(x0))
  dots <- plot2Drglbis("segments3D", x0 = x0, y0 = y0, x1 = x1, y1 = y1,
    dz = dz, colkey = colkey, colvar = colvar, col = col,
    NAcol = NAcol, breaks = breaks, clim = clim, clab = clab, add = add,
    namesextra = namesextra, ...)
  finishplotrgl(dots, namesextra)
}
