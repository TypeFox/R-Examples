 mds <- function(delta, ndim = 2, type = c("ratio", "interval", "ordinal","mspline"), 
                         weightmat = NULL, init = NULL, ties = "primary",  verbose = FALSE, 
                         relax = FALSE, modulus = 1, itmax = 1000, eps = 1e-6, 
                         spline.degree = 2, spline.intKnots = 2) {
   smacofSym(delta=delta, ndim = ndim, type = type, 
             weightmat = weightmat, init = init, ties = ties,  verbose = verbose, 
             relax = relax, modulus = modulus, itmax = itmax, eps = eps, 
             spline.degree = spline.degree, spline.intKnots = spline.intKnots)
 }
   
#mds <- smacofSym