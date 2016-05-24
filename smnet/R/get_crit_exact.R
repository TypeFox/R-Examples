get_crit_exact<-function(rhoArgs, X, XTX, P, response, cholFactor, n, 
                          np = nrow(XTX), Xw, Xy, n.sm, identifyBit, crit){    

  for(j in 1:n.sm) P[[j]]<-P[[j]]*exp(rhoArgs[j])
  Psum     <- Reduce("+", P)
  info     <- XTX + Psum + identifyBit
  U        <- try(update.spam.chol.NgPeyton(cholFactor, info), silent = T)
  if(class(U) == "try-error"){
    out <- 10^20
  } else {
    beta_hat <- backsolve.spam(U, forwardsolve.spam(U, Xy))  
    left1    <- forwardsolve.spam(U, t(X))
    ED1      <- sum(left1*left1)
    ED       <- ED1
    fit      <- X %*% beta_hat  
    if(crit == "AICC") out <- log(sum((response - fit)^2)) + (2*(ED+1)/(n-ED-2))
    if(crit == "AIC")  out <- log(sum((response - fit)^2)) + 2*ED/n
    if(crit == "GCV")  out <- sum((response - fit)^2)/(1-(ED/n))^2
  }
  out
}
