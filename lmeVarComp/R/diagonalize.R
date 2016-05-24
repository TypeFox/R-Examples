diagonalize <-
function(x, tol = 1e-10)
{
  p <- nrow(x)
  k <- ncol(x) %/% nrow(x)
  
  result <- .C("diagonalize", A = as.double(x), 
    p = as.integer(p), k = as.integer(k), tol = as.double(tol), 
    D = double(p * k), fail = integer(1L),
    DUP = TRUE)
  
  list(U = matrix(result$A, p, p), D = matrix(result$D, p, k))
}
