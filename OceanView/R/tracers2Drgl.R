
## =============================================================================
## =============================================================================
## Tracers in 2D using rgl
## =============================================================================
## =============================================================================

tracers2Drgl <- function(x, y, colvar = NULL, ..., 
  col = NULL, NAcol = "white", breaks = NULL,
  colkey = FALSE, clim = NULL, clab = NULL) {
   z <- rep (1.001, length.out = length(x))
   tracers3Drgl(x, y, z, colvar = colvar, 
     col = col, NAcol = NAcol, breaks = breaks, colkey = colkey,
     clim = clim, clab = clab, ...)
}

