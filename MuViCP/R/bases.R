normalise <- function(x) {
  if (is.matrix(x)) {
    lengths <- sqrt(colSums(x ^ 2, na.rm = TRUE))
    sweep(x, 2, lengths, "/")
  } else {
    x / sqrt(sum(x ^ 2))    
  }
}

orthonormalise <- function(x) {
  x <- normalise(x) # to be conservative
  
  if (ncol(x) > 1) {
    for (j in seq_len(ncol(x))) {
      for (i in seq_len(j - 1)) {
        x[, j] <- x[, j] - crossprod(x[, j], x[, i]) * x[, i]
      }
    }
  }
  
  normalise(x)
}

basis_random <- function (n, d = 2) 
{
    mvn <- matrix(rnorm(n * d), ncol = d)
    orthonormalise(mvn)
}
