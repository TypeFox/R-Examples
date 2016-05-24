perm <- 
function(I, n, i) {
  n <- n - 1
  p <- c()
  while (length(I) > 0) {
    f <- factorial(length(I)-1)
    j <- floor(n / f)
    x <- I[j+1]
    n <- n %% f
    p <- c(p, x)
    I <- I[-j-1]
  }
  return(c(p, i))
}