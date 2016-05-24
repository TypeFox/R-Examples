#'# mExp --- calculate the exponent of a given matrix
#'# INPUT = a matrix (x) and the exponent (n)
#'# OUTPUT = the resulting exponentiated matrix
#'# 
#'# Alberto Monteiro (https://stat.ethz.ch/pipermail/
#'# r-help/2007-May/131330.html)
#'# ___________________________________________________

mExp <- function(x='matrix', n=2){
  if (n == 1) return(x)
  result <- diag(1, ncol(x))
  while (n > 0) {
    if (n %% 2 != 0) {
      result <- result %*% x
      n <- n - 1
    }
    x <- x %*% x
    n <- n / 2
  }

  rownames(result) <- colnames(result)

  return(result)

}
