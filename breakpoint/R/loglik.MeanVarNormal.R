loglik.MeanVarNormal <-
function(a, b, data, h){
  if((b-a) <= h) {
    ll <- -100000000
    
  } else {
    x <- data[a:(b - 1), 1]
    mu <- mean(x)    
    v <- var(x)
    
    ll <- try(- {1/2 * length(x)*log(2*pi*v)} - {1/2* 1/v * sum((x - mean(x))^2)}, silent=T)
    
    if(mode(ll) == "character"){ 
      ll <- -10000000
    } else ll <- ll
  }
  return(list(ll = ll))
}
