

fitclogitboostpredict <- function(fit, M, groupid){
  addfunction <- fit$func
  addfunction.var <- fit$index
  theta <- fit$theta
  rho <- fit$rho
  bigf <- function(x){
    temp <- 0
    for (m in seq(1, length(theta))){
      temp <- temp + rho * theta[[m]] * predict(addfunction[[m]], x[addfunction.var[[m]]])$y
    }
    return(temp)
  }
  
  fx <- apply(M, 1, bigf)  
  expprob <- exp(fx)
  expprob.sum <- aggregate(expprob,by=list(groupid),FUN=sum)
  
  probest <- rep(0, length(fx))
  for (i in seq(1, length(fx))){
    probest[i] <- expprob[i] / expprob.sum[expprob.sum[, 1] == groupid[i], 2]
  }
  
  return(list(prob = probest, utility = fx))
}