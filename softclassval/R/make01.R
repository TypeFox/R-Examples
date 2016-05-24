##' @include softclassval.R
## helper for hard "and" operator: sets all values NA that are not within +- tol from 0 or 1.
.make01 <- function (x, tol = 1e-6){
  tmp <- rep (NA_real_, length (x))
  tmp [x >    -tol & x <     tol] <- 0
  tmp [x > 1 - tol & x < 1 + tol] <- 1
  
  attributes (tmp) <- attributes (x)

  tmp
}
.test (.make01) <- function (){
  checkIdentical (.make01 (v), c( a = 0, b = NA, c = NA, d = 1, e = NA))
  checkIdentical (attributes (.make01 (m)), attributes (m))
}
