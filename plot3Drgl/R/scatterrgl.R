# ==============================================================================
# the 2D scatterplot (lines, points) function, using rgl
# ==============================================================================
lines2Drgl <- function(x, y, ...) {
    dot <- list(...)
    if (is.null(dot$type))
        dot$type <- "l"
    do.call("scatter2Drgl", c(alist(x, y), dot))
}

points2Drgl <- function(x, y, ...) {
    dot <- list(...)
    if (is.null(dot$type))
        dot$type <- "p"
    do.call("scatter2Drgl", c(alist(x, y), dot))
}


scatter2Drgl <- function(x, y, colvar = NULL, ...,
              col = NULL, NAcol = "white", breaks = NULL, colkey = NULL,
              clim = NULL, clab = NULL, CI = NULL, dz = 0.1, add = FALSE)  {
# ------------------------------------------------------------------------------

  namesextra <- c("type", "pch", "cex", "lwd", "lty", "CI")
  type <- list(...)$type
  if (is.null(type))
    type <- "p"

  dots <- plot2Drgl("scatter3D", x, y,
     colvar, col, NAcol, breaks = breaks,
     clim, add, clab = clab, z = rep(1 + dz, length.out = length(x)),
     namesextra = namesextra, CI = CI, colkey = colkey, ...)
  if (type == "h") {
    ymin <- max(c(min(y), 0))
    segments2Drgl(x, rep(ymin, length.out = length(x)), x, y,
     colvar = colvar, col = col, NAcol = NAcol, breaks = breaks,
     clim = clim, add = TRUE, ...)
  }
  finishplotrgl(dots, namesextra = namesextra)
}

# ==============================================================================
# the 2D text function, using rgl
# ==============================================================================

text2Drgl <- function(x, y, labels, colvar = NULL, ...,
             col = NULL, NAcol = "white", breaks = NULL, colkey = NULL,
             clim = NULL, clab = NULL, dz = 0.1, add = FALSE)  {
# ------------------------------------------------------------------------------

  namesextra <- c("labels", "cex", "font")
  dots <- plot2Drgl("text3D", x, y, z = rep(1 + dz, length.out = length(x)),
     colvar, col, NAcol, breaks = breaks, clim, add, clab = clab,
     namesextra = namesextra, labels = labels, colkey = colkey, ...)
  finishplotrgl(dots, namesextra = namesextra)
}
