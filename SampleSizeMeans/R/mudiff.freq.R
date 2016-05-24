`mudiff.freq` <-
function(len,lambda1,lambda2,level=0.95,equal=TRUE)
{
  z <- qnorm((1+level)/2)
  
  if (!equal)
  {
    n1 <- 4*z^2/len^2*(1/sqrt(lambda1*lambda2)+1/lambda1)
    n2 <- sqrt(lambda1/lambda2)*n1
    n1 <- floor(n1) + rep(c(0,1),c(2,2))
    n2 <- floor(n2) + rep(c(0, 1), 2)
    v <- 1/(lambda1*n1) + 1/(lambda2*n2)
    w <- which(2*z*sqrt(v) <= len)
    n1 <- n1[w]
    n2 <- n2[w]
    tot.n <- n1 + n2
    w <- which(tot.n==min(tot.n))
    n <- c(n1[w[1]], n2[w[1]])
  }
  else
  {
    n <- rep(ceiling(4*z^2/len^2*(1/lambda1+1/lambda2)), 2)
  }

  n
}

