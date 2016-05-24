
#' @importFrom graphics lines polygon

curveseg <- function(x0, x1, y0, y1, width = 1, nsteps = 50,
                     colorstyle, col = "#ffcc0066", grad = NULL, lty = 1,
                     curvestyle = c("sin", "line")) {

  curvestyle <- match.arg(curvestyle)

  w <- width

  if (colorstyle == "gradient") {
    grad <- color_ramp_palette_alpha(grad)(nsteps)

  } else {
    grad <- rep(col, nsteps)
  }

  if (curvestyle == "sin" ) {
    xx <- seq(-pi/2, pi/2, length.out = nsteps)
    yy <- y0 + (y1 - y0) * (sin(xx) + 1 ) / 2
    xx <- seq(x0, x1, length.out = nsteps)
  }

  if (curvestyle == "line" ) {
    xx <- seq(x0, x1, length.out = nsteps)
    yy <- seq(y0, y1, length.out = nsteps)
  }

  for(i in 1:(nsteps-1) ) {
    polygon(
      c(xx[i], xx[i + 1], xx[i + 1], xx[i] ),
      c(yy[i], yy[i + 1], yy[i + 1] + w, yy[i] + w ),
      col = grad[i], border = grad[i]
    )

    lines(c(xx[i], xx[i + 1]), c(yy[i], yy[i + 1]), lty = lty)
    lines(c(xx[i], xx[i + 1]), c(yy[i] + w, yy[i + 1] + w ), lty = lty)
  }
}
