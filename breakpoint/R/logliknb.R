logliknb <-
function(a, b, data, r, h){  
  if((b-a) <= h) {
    ll <- -100000000
  } else {
    x <- data[a:(b - 1), 1]
    mu <- mean(x)    
    ll <- try(sum(lgamma(x + r)) - sum(lgamma(x + 1)) - length(x) * lgamma(r) + log(mu/r) * sum(x) - log(1 + (mu / r)) * sum(x + r), silent=T)
    if(mode(ll) == "character"){ 
      ll <- -10000000
    } else {
      ll <- ll
    }
  }
  return(list(ll = ll))
}
