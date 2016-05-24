## sanity checks for dissimilatiries (>=0)

checkdiss <- function(diss) {
  if (any(sapply(diss, function(d0) any(d0 < 0, na.rm = TRUE)))) stop("Dissimilarities should be non-negative!")
}