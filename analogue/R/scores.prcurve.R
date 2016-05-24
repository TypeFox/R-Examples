## `scores` extractor function for prcurve class
`scores.prcurve` <- function(x, display = c("curve","dimensions"), ...) {
  display <- match.arg(display)
  ## return position along the curve?
  if (isTRUE(all.equal(display, "curve"))) {
    ret <- matrix(x$lambda, ncol = 1)
    rownames(ret) <- names(x$lambda)
    colnames(ret) <- "PrC"
  }
  ## return coordinates of curve in each dimension
  if (isTRUE(all.equal(display, "dimensions"))) {
    ret <- x$s
  }
  ret
}
