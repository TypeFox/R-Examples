sinv <-
function(x)
{
  n <- ncol(x)
  
  result <- .C("sinv", x = as.double(x), n = as.integer(n),
    DUP = TRUE)
  
  matrix(result$x, n, n)
}
