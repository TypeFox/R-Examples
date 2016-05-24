loglikzinb <-
function(a, b, data, r, h){    
  x <- data[a : (b - 1), 1]
  sumzero <- sum( x== 0)
  mu <- mean(x)
  
  if((b-a) <= h) {
    ll<- -10000000
  
  } else if (sumzero == 0) {
    ll <- try(sum(lgamma(x + r)) - sum(lgamma(x + 1)) - length(x) * lgamma(r) + log(mu / r) * sum(x) - log(1 + (mu / r)) * sum(x + r), silent = T)
    if(mode(ll) == "character") {
      ll <- -10000000
    } else {
      ll<-ll
    }
  } else {
    p0 <- sumzero / length(x)
    ll1 <- try(sumzero * log(p0 + (1 - p0) * (r / (r + mu)) ^ r), silent = T) 
      if(mode(ll1) == "character") {
        ll1 <- -10000000
    } else {
      ll1 <- ll1
    }
    newv <- x[x > 0]
    Ns <- length(newv)
    
    ll2 <- try(Ns * log(1 - p0) + sum(lgamma(newv + r)) - sum(lgamma(newv + 1)) - Ns * lgamma(r) + log(mu / r) * sum(newv) - log(1 + (mu / r)) * sum(newv + r), silent = T)
    if(mode(ll2) == "character") {
      ll2 <- -10000000
    } else {
      ll2 <- ll2
    }
    ll <- ll1 + ll2 
  }   
  return(list(ll = ll))
}
