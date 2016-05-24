`mu.modwoc` <-
function(len, alpha, beta, n0, level = 0.95, worst.level = 0.95)
{
  l.prior <- 2*sqrt(beta/alpha/n0)*qt((1+level)/2,2*alpha)
  
  if (l.prior <= len) {
    0 # prior knowledge is sufficient
  }
  else {
    step <- 2^6
    found.upper.bound <- FALSE
    found.lower.bound <- FALSE
    direction <- 1
    n <- 8*beta/len^2*qnorm((1+level)/2)^2/qchisq(1-worst.level,2*alpha)
    n <- ceiling(max(1,n-n0))

    while(step >= 2) {
      step <- ifelse(found.upper.bound & found.lower.bound,step/2,step*2)
      n <- n + direction * step

      if(n <= 0) {
        found.lower.bound <- TRUE
        n <- 1
      }

      a <- qf(worst.level, n, 2 * alpha)
      denom <- (1 + n/2/alpha * a)*8*beta
      num <- (n+n0) * len^2 *(2*alpha+n)
      b <- pf(num/denom,1,n+2*alpha)
      if(b >= level) {
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

