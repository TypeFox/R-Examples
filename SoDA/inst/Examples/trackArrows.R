arrowFrom <- function(u, fraction) {
      n = length(u)
      if(n < 2)
          numeric(0)
      else
          u[-1]*(1-fraction) + u[-n]*fraction
  }

trackArrows <- function(x, y,
                        fraction, head, nArrowLengths = 5, ...) {
    x0 = arrowFrom(x, fraction);  y0 = arrowFrom(y, fraction)
    x1 = x[-1];   y1 = y[-1]
    ## compute the average line length
    delta = sqrt(mean((x1-x0)^2 + (y1-y0)^2, na.rm = TRUE))
    ## and convert it to inches for arrows()
    delta = delta * (par("pin")[1]/diff(range(x, na.rm = TRUE)))
    arrows(x0, y0, x1, y1, head * delta, ...)
}
