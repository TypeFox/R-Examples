`mu.alc` <-
function(len,alpha,beta,n0,level=0.95)
{
  l.prior <- 2*sqrt(beta/alpha/n0)*qt((1+level)/2,2*alpha)
  
  if (l.prior <= len)
  {
    0 # prior knowledge is sufficient
  }
  else
  {
    step <- 2^3
    found.upper.bound <- FALSE
    found.lower.bound <- FALSE
    direction <- 1
    bn <- 2*qnorm((1+level)/2)*sqrt(beta)/len*exp(lgamma(alpha-.5)-lgamma(alpha))
    n <- ceiling(max(1,bn^2-n0))

    while(step >= 2) {
      step <- ifelse(found.upper.bound & found.lower.bound, step/2, step*2)
      n <- n + direction * step

      if(n <= 0) {
        found.lower.bound <- TRUE
        n <- 1
      }

      under.root <- 2*beta/(2*alpha+n)/(n+n0)
      E <- exp(lgamma(n/2+alpha)-lgamma((n-1)/2+alpha))
      E <- E*exp(lgamma(alpha-1/2)-lgamma(alpha))
      avg.len <- 2*qt((1+level)/2,n+2*alpha)*sqrt(under.root)*E

      if(avg.len<=len) {
        found.upper.bound <- TRUE
        direction <- -1
      }
      else {
        found.lower.bound <- TRUE
        direction <- 1
      }

    }
    n[direction == 1] <- n + 1
    n
  }
}

