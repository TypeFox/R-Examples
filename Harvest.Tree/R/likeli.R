likeli <- function(n, m){ 
  if (m < n & m != 0)
    Like <- m*log(m/n) + (n-m)*log(1-m/n)
  else Like <- 0
  return(Like)
}