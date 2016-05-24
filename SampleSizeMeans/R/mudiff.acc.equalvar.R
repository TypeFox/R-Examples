`mudiff.acc.equalvar` <-
function(len,alpha,beta,n01,n02,level=0.95,equal=TRUE)
{
  l.prior <- 2*qt((1+level)/2,2*alpha)*sqrt(beta/alpha*(1/n01+1/n02))
  if (l.prior <= len){
    0 # prior knowledge is sufficient
  }
  else {
    t2 <- qt((1+level)/2,2*alpha)^2

    if (!equal){
      m <- t2*8*beta/alpha/len^2
      N <- ceiling(m)-c(n01,n02)

      # if one of the prior sample sizes is sufficient...
      d <- 4*beta/alpha/len/len*t2
      if (N[2] < 0)
      {
        N[1] <- ceiling(1/(1/d-1/n02)-n01)
        N <- pmax(0,N)
      }
      if (N[1] < 0)
      {
        N[2] <- ceiling(1/(1/d-1/n01)-n02)
        N <- pmax(0,N)
      }
    }
    else {
      A <- alpha*len*len/4
      B <- A*(n01+n02)-2*beta*t2
      C <- n01*n02*A-beta*t2*(n01+n02)

      N <- rep(quadsol(A,B,C),2)
    }

    N
  }
}

