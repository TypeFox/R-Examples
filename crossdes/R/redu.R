redu <- function(a, b) {
  
  la <- length(a)
  lb <- length(b)
  if (la < lb || lb < 2) {
    stop("Please check your polynomials.")
  }
  
  b1 <- b/b[1]  # divide by coefficient of highest power
  
  dummy <- rep(0, la)  # partial result to be subtracted
  
  for (i in 1:(la - lb + 1)) {
    dummy <- numeric(la)
    dummy[i:(i + lb - 1)] <- a[i] * b1
    a <- a - dummy
  }
  a[(la - lb + 2):la]
  
} 
