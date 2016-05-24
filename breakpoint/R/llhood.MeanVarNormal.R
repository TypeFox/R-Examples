llhood.MeanVarNormal <-
function(y, data, h){
  llihood <- c()
  for(j in 1:(length(y) - 1)){
    llihood[j] <- loglik.MeanVarNormal(y[j], y[j + 1], data, h)$ll  
  }
  loglikelihood <- sum(llihood)
}
