
#' @importFrom grDevices rgb

color_ramp_palette_alpha <- function (colors, ...) {
  ramp <- color_ramp_alpha(colors, ...)
  function(n) {
    x <- ramp(seq.int(0, 1, length.out = n))
    rgb(x[, 1L], x[, 2L], x[, 3L], x[ , 4L], maxColorValue = 255)
  }
}

#' @importFrom grDevices col2rgb
#' @importFrom stats approxfun splinefun

color_ramp_alpha <- function (colors, bias = 1,
                              interpolate = c("linear", "spline")) {

  if (bias <= 0) stop("'bias' must be positive")
  interpolate <- match.arg(interpolate)
  colors <- t(col2rgb(colors, alpha = TRUE) / 255)

  interpolate <- switch(
    interpolate,
    linear = stats::approxfun,
    spline = stats::splinefun
  )

  if ((nc <- nrow(colors)) == 1L) {
    colors <- colors[c(1L, 1L), ]
    nc <- 2L
  }

  x <- seq.int(0, 1, length.out = nc) ^ bias

  palette <- c(
    interpolate(x, colors[, 1L]),
    interpolate(x, colors[, 2L]),
    interpolate(x, colors[, 3L]),
    interpolate(x, colors[, 4L])
  )

  roundcolor <- function(rgb) pmax(pmin(rgb, 1), 0)

  function(x) 255 * roundcolor(cbind(
    palette[[1L]](x),
    palette[[2L]](x),
    palette[[3L]](x),
    palette[[4L]](x)
  ))
}
