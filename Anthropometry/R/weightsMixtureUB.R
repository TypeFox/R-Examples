weightsMixtureUB <- function(orness,numVar){
  if(orness == .5){
    w = rep(1/numVar,numVar)
   }else{
     lambda = 0.5
     prob0 = 3/2 - 2 * orness
     w = lambda * dbinom(0:(numVar-1), size = numVar-1, prob = prob0) + 
       (1 - lambda) * rep(1/numVar,numVar)
    }
  w
}

