`mu.mblalc` <-
function(len,alpha,beta,level=0.95)
{
  n <- ceiling(4*beta*(qnorm((1+level)/2)/len*exp(lgamma(alpha-.5)-lgamma(alpha)))^2)
  step <- 2^3

  found.upper.bound <- FALSE
  found.lower.bound <- FALSE
  direction <- +1

  while(step >= 2)
  {
          step <- ifelse(found.upper.bound & found.lower.bound,
                      step/2,step*2)
          n <- n + direction * step
          if(n <= 2) {
            found.lower.bound <- TRUE
            n <- 2
          }

          fn <- 2*qt((1+level)/2,n-1)*sqrt(2*beta/n/(n-1))*
             exp(lgamma(alpha-.5)-lgamma(alpha))*exp(lgamma(n/2)-lgamma((n-1)/2))

          if(fn<=len) {
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

