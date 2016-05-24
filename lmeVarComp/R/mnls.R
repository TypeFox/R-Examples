mnls <-
function(x, y, rcond = 1e-10)
{
  m <- nrow(x)
  n <- ncol(x)
  nrhs <- NCOL(y)
  
  result <- .C("mnls", x = as.double(x), y = as.double(y), 
    beta = double(n * nrhs), m = as.integer(m), n = as.integer(n),
    nrhs = as.integer(nrhs), rcond = as.double(rcond), 
    rank = integer(1L), 
    DUP = TRUE)
  
  beta <- matrix(result$beta, n, nrhs)
  attr(beta, "rank") <- result$rank
  beta
}
