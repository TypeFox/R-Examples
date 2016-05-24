cmd <- function(m1, m2) {
  if (!is.matrix(m1) || !is.matrix(m2))
    stop("requires two matrices to compute the cmd")

  if (any(dim(m1) - dim(m2) != 0))
    stop("the two matrices must have the same dimensions")
  
  1 - sum(diag(m1 %*% m2)) /(norm(m1, type="F")*norm(m2, type="F"))
}
