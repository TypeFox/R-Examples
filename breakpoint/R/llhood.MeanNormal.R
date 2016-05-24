llhood.MeanNormal <-
function(y, data, v, h){
  llihood <- c()
  for(j in 1:(length(y) - 1)){
    llihood[j] <- loglik.MeanNormal(y[j], y[j + 1], data, h, v)$ll  
  }
  loglikelihood <- sum(llihood)
}
