llhoodnb <-
function(y, data, r, h){
  llihood <- c()
    for(j in 1:(length(y) - 1)){
      llihood[j] <- logliknb(y[j], y[j + 1], data, r, h)$ll  
    }
  loglikelihood <- sum(llihood)
}
