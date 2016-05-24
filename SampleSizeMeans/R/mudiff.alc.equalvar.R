`mudiff.alc.equalvar` <-
function(len,alpha,beta,n01,n02,level=0.95,equal=TRUE)
{
  l.prior <- 2*qt((1+level)/2,2*alpha)*sqrt(beta/alpha*(1/n01+1/n02))
  if (l.prior <= len) {
    0 # prior knowledge is sufficient
  }
  else {
    step <- 2^3
    found.upper.bound <- FALSE
    found.lower.bound <- FALSE
    direction <- 1
    bn <- 2*qnorm((1+level)/2)*sqrt(2*beta)/len*exp(lgamma(alpha-.5)-lgamma(alpha))
    n1 <- ceiling(max(0,bn^2-n01))

    while(step >= 2) {
      step <- ifelse(found.upper.bound & found.lower.bound, step/2, step*2)
      n1 <- n1 + direction * step
      if(n1 <= 0) {
        found.lower.bound <- TRUE
        n1 <- 0
      }
      n2 <- ifelse(!equal,n1+n01-n02,n1)
      n2 <- max(0,n2)
      ndot <- n1+n2
      under.root <- 2*beta*(1/(n1+n01)+1/(n2+n02))/(2*alpha+ndot)
      E <- exp(lgamma(ndot/2+alpha)-lgamma((ndot-1)/2+alpha))
      E <- E*exp(lgamma(alpha-1/2)-lgamma(alpha))
      avg.len <- 2*qt((1+level)/2,ndot+2*alpha)*sqrt(under.root)*E

      if(avg.len<=len) {
        found.upper.bound <- TRUE
        direction <- -1
      }
      else {
        found.lower.bound <- TRUE
        direction <- 1
      }
    }
    n1[direction == 1] <- n1 + 1
    n2 <- ifelse(!equal,n1+n01-n02,n1)
    n2 <- max(0,n2)

    c(n1,n2)
  }
}

