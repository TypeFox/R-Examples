critAgcv <- function(k,valpr,tPADmdemiY,DdemiPA,n,ddlmini,ddlmaxi){
  ImLambdak <- rep(0,n)
  if (ddlmini>=1) {
    valpr0 <- valpr[-(1:ddlmini)]
    ImLambdak[-(1:ddlmini)] <- (1-valpr0)^k
  } else {
    ImLambdak <- (1-valpr)^k
  }   
  sigma2 <- sum((DdemiPA%*%(ImLambdak*tPADmdemiY))^2)/n
  logsigma2 <- log(sigma2)
  traceSk <- sum(1-ImLambdak)
  if (traceSk>ddlmaxi)  return(.Machine$double.xmax)
  if ((sigma2<1e-10)|(traceSk>((1-1e-10)*n))) return(.Machine$double.xmax)
  return(logsigma2 - 2*log(1-traceSk/n))
}

critAaicc <- function(k,valpr,tPADmdemiY,DdemiPA,n,ddlmini,ddlmaxi){
  ImLambdak <- rep(0,n)
  if (ddlmini>=1) {
    valpr0 <- valpr[-(1:ddlmini)]
    ImLambdak[-(1:ddlmini)] <- (1-valpr0)^k
  } else {
    ImLambdak <- (1-valpr)^k
  }   
  sigma2 <- sum((DdemiPA%*%(ImLambdak*tPADmdemiY))^2)/n
  logsigma2 <- log(sigma2)
  traceSk <- sum(1-ImLambdak)
  if (traceSk>ddlmaxi)  return(.Machine$double.xmax)
  if ((sigma2<1e-10)|(traceSk>((1-1e-10)*n))) return(.Machine$double.xmax)
  return(logsigma2 + 1+ 2*(traceSk+1)/(n-traceSk-2))
}

critAaic <- function(k,valpr,tPADmdemiY,DdemiPA,n,ddlmini,ddlmaxi){
  ImLambdak <- rep(0,n)
  if (ddlmini>=1) {
    valpr0 <- valpr[-(1:ddlmini)]
    ImLambdak[-(1:ddlmini)] <- (1-valpr0)^k
  } else {
    ImLambdak <- (1-valpr)^k
  }   
  sigma2 <- sum((DdemiPA%*%(ImLambdak*tPADmdemiY))^2)/n
  logsigma2 <- log(sigma2)
  traceSk <- sum(1-ImLambdak)
  if (traceSk>ddlmaxi)  return(.Machine$double.xmax)
  if ((sigma2<1e-10)|(traceSk>((1-1e-10)*n))) return(.Machine$double.xmax)
  return(logsigma2 + 2*(traceSk)/(n))
}
critAbic <- function(k,valpr,tPADmdemiY,DdemiPA,n,ddlmini,ddlmaxi){
  ImLambdak <- rep(0,n)
  if (ddlmini>=1) {
    valpr0 <- valpr[-(1:ddlmini)]
    ImLambdak[-(1:ddlmini)] <- (1-valpr0)^k
  } else {
    ImLambdak <- (1-valpr)^k
  }   
  sigma2 <- sum((DdemiPA%*%(ImLambdak*tPADmdemiY))^2)/n
  logsigma2 <- log(sigma2)
  traceSk <- sum(1-ImLambdak)
  if (traceSk>ddlmaxi)  return(.Machine$double.xmax)
  if ((sigma2<1e-10)|(traceSk>((1-1e-10)*n))) return(.Machine$double.xmax)
  return(logsigma2 + log(n)*(traceSk)/(n))
}
critAgmdl <- function(k,valpr,tPADmdemiY,DdemiPA,n,ddlmini,ddlmaxi,Y){
  ImLambdak <- rep(0,n)
  if (ddlmini>=1) {
    valpr0 <- valpr[-(1:ddlmini)]
    ImLambdak[-(1:ddlmini)] <- (1-valpr0)^k
  } else {
    ImLambdak <- (1-valpr)^k
  }   
  sigma2 <- sum((DdemiPA%*%(ImLambdak*tPADmdemiY))^2)/n
  logsigma2 <- log(sigma2)
  traceSk <- sum(1-ImLambdak)
  Sbul <- n*sigma2/(n-traceSk)
  if (traceSk>ddlmaxi)  return(.Machine$double.xmax)
  if ((sigma2<1e-10)|(traceSk>((1-1e-10)*n))) return(.Machine$double.xmax)
  return(log(Sbul)+traceSk/n*log((sum(Y^2)-n*sigma2)/(traceSk*Sbul)))
}

critS1gcv <- function(k,valpr,tUy2,n,ddlmini,ddlmaxi){
  valpr0 <- valpr[-(1:ddlmini)]
  ImLambdak <- rep(0,n)
  ImLambdak[-(1:ddlmini)] <- (1-valpr0)^k
  logsigma2 <- log(sum((ImLambdak^2)*(tUy2)))-log(n)  
  traceSk <- sum(1-ImLambdak)
  if (traceSk>ddlmaxi)  return(.Machine$double.xmax)
  if ((logsigma2<log(1e-10))|(traceSk>((1-1e-10)*n))) return(.Machine$double.xmax)
  return(logsigma2 - 2*log(1-traceSk/n))
}

critS1aicc <- function(k,valpr,tUy2,n,ddlmini,ddlmaxi){
  valpr0 <- valpr[-(1:ddlmini)]
  ImLambdak <- rep(0,n)
  ImLambdak[-(1:ddlmini)] <- (1-valpr0)^k
  logsigma2 <- log(sum((ImLambdak^2)*(tUy2)))-log(n)  
  traceSk <- sum(1-ImLambdak)
  if (traceSk>ddlmaxi)  return(.Machine$double.xmax)
  if ((logsigma2<log(1e-10))|(traceSk>((1-1e-10)*n))) return(.Machine$double.xmax)
  return(logsigma2 + 1+ 2*(traceSk+1)/(n-traceSk-2))
}

critS1aic <- function(k,valpr,tUy2,n,ddlmini,ddlmaxi){
  valpr0 <- valpr[-(1:ddlmini)]
  ImLambdak <- rep(0,n)
  ImLambdak[-(1:ddlmini)] <- (1-valpr0)^k
  logsigma2 <- log(sum((ImLambdak^2)*(tUy2)))-log(n)  
  traceSk <- sum(1-ImLambdak)
  if (traceSk>ddlmaxi)  return(.Machine$double.xmax)
  if ((logsigma2<log(1e-10))|(traceSk>((1-1e-10)*n))) return(.Machine$double.xmax)
  return(logsigma2 + 2*(traceSk)/(n))
}
critS1bic <- function(k,valpr,tUy2,n,ddlmini,ddlmaxi){
  valpr0 <- valpr[-(1:ddlmini)]
  ImLambdak <- rep(0,n)
  ImLambdak[-(1:ddlmini)] <- (1-valpr0)^k
  logsigma2 <- log(sum((ImLambdak^2)*(tUy2)))-log(n)  
  traceSk <- sum(1-ImLambdak)
  if (traceSk>ddlmaxi)  return(.Machine$double.xmax)
  if ((logsigma2<log(1e-10))|(traceSk>((1-1e-10)*n))) return(.Machine$double.xmax)
  return(logsigma2 + log(n)*(traceSk)/(n))
}
critS1gmdl <- function(k,valpr,tUy2,n,ddlmini,ddlmaxi,Y){
  valpr0 <- valpr[-(1:ddlmini)]
  ImLambdak <- rep(0,n)
  ImLambdak[-(1:ddlmini)] <- (1-valpr0)^k
  sigma2 <- sum((ImLambdak^2)*(tUy2))/n
  logsigma2 <- log(sigma2)
  traceSk <- sum(1-ImLambdak)
    Sbul <- n*sigma2/(n-traceSk)
  if (traceSk>ddlmaxi)  return(.Machine$double.xmax)
  if ((logsigma2<log(1e-10))|(traceSk>((1-1e-10)*n))) return(.Machine$double.xmax)
    return(log(Sbul)+traceSk/n*log((sum(Y^2)-n*sigma2)/(traceSk*Sbul)))
}

 
 
