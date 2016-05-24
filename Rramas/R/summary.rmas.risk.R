summary.rmas.risk <-
function(object, q=c(0.025,0.975),...){
  cosa <- object
  cosa.boot <- NULL 

    # extract bootstraped probabilities
  for(i in 1:length(cosa$cf.boot)){
    cosa.boot <- cbind(cosa.boot,cosa$cf.boot[[i]][,2])
  }
    # paste together Threshold, Probanbility (cosa$cf.obs) and 95% C.I.
  tabla <- cbind(cosa$cf.obs,t(apply(cosa.boot,1,quantile, q)))
  class(tabla) <- c("summary.rmas.risk", class(tabla))
    
  plot(tabla, main= cosa$main)
  return(tabla)
}

