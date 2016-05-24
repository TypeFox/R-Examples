indscal <- function(delta, ndim = 2, type = c("ratio", "interval", "ordinal", "mspline"), 
                   weightmat = NULL, init = "torgerson", ties = "primary", 
                   verbose = FALSE, modulus = 1, itmax = 1000, eps = 1e-6,
                   spline.degree = 2, spline.intKnots = 2) {
  
  smacofIndDiff(delta = delta, ndim = ndim, type = type, constraint = "indscal",
                weightmat = weightmat, init = init, ties = ties, verbose = verbose, modulus = modulus, 
                itmax = itmax, eps = eps, spline.degree = spline.degree, spline.intKnots = spline.intKnots)
}