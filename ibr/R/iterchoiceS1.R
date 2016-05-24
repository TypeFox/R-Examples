iterchoiceS1 <- function(n,mini,maxi,tUy,eigenvaluesS1,ddlmini,ddlmaxi,y,criterion,fraction){
  res <- list(minimum=0,objective=.Machine$double.xmax)
  tUy2 <- tUy^2
  if (criterion=="gmdl") {
    objectif <- .Machine$double.xmax
    fraction <- sort(unique(c(fraction[(fraction<maxi)&(fraction>mini)],mini,maxi)))
    dep <- length(fraction)
    repeat {
      mini <- fraction[dep-1]
      objectif <- critS1gmdl(mini,valpr=eigenvaluesS1,tUy2=tUy2,n=n,ddlmini=ddlmini,ddlmaxi=ddlmaxi,Y=y)
      if (objectif<.Machine$double.xmax) {
        bb <- fraction[dep]
        break 
      }
      dep <- dep-1
      if (dep==1) stop(paste("decrease Kmax below",fraction[2],"or increase dfmaxi component of control.par list"))
    }
    phi <- (sqrt(5) - 1)/2
    repeat {
      x1 <- mini + (1-phi)*(bb-mini)
      objectif <- critS1gmdl(x1,valpr=eigenvaluesS1,tUy2=tUy2,n=n,ddlmini=ddlmini,ddlmaxi=ddlmaxi,Y=y)
      if (objectif<.Machine$double.xmax) break else {
        bb <- x1
      }
    }
    fraction <- c(fraction[1:(dep-1)],bb)
    for (i in 1:(length(fraction)-1)) {
      res1 <- optimize(critS1gmdl,lower=fraction[i],upper=fraction[i+1],tol=0.5,valpr=eigenvaluesS1,tUy2=tUy2,n=n,ddlmini=ddlmini,ddlmaxi=ddlmaxi,Y=y)
      if (res1$objective<res$objective) res <- res1
    }
  } else {
    fcriterion <- switch(criterion,aic=critS1aic,aicc=critS1aicc,gcv=critS1gcv,bic=critS1bic)
    objectif <- .Machine$double.xmax
    fraction <- sort(unique(c(fraction[(fraction<maxi)&(fraction>mini)],mini,maxi)))
    dep <- length(fraction)
    repeat {
      mini <- fraction[dep-1]
      objectif <- fcriterion(mini,valpr=eigenvaluesS1,tUy2=tUy2,n=n,ddlmini=ddlmini,ddlmaxi=ddlmaxi)
      if (objectif<.Machine$double.xmax) {
        bb <- fraction[dep]
        break 
      }
      dep <- dep-1
      if (dep==1) stop(paste("decrease Kmax below",fraction[2],"or increase dfmaxi component of control.par list"))
    }
    phi <- (sqrt(5) - 1)/2
    repeat {
      x1 <- mini + (1-phi)*(bb-mini)
      objectif <- fcriterion(x1,valpr=eigenvaluesS1,tUy2=tUy2,n=n,ddlmini=ddlmini,ddlmaxi=ddlmaxi)
      if (objectif<.Machine$double.xmax) break else {
        bb <- x1
      }
    }
    fraction <- c(fraction[1:(dep-1)],bb)
    for (i in 1:(length(fraction)-1)) {
      res1 <- optimize(fcriterion,lower=fraction[i],upper=fraction[i+1],tol=0.5,valpr=eigenvaluesS1,tUy2=tUy2,n=n,ddlmini=ddlmini,ddlmaxi=ddlmaxi)
      if (res1$objective<res$objective) res <- res1
    }     
  }
  return(list(iter=floor(res$minimum),objective=res$objective))
}    

iterchoiceS1e <- function(y,K,tUy,eigenvaluesS1,ddlmini,ddlmaxi){
  n <- length(y)
  iter <- 1
  tUy2 <- tUy^2
  GCV <- rep(.Machine$double.xmax,length(K))
  AICc <- rep(.Machine$double.xmax,length(K))
  AIC <- rep(.Machine$double.xmax,length(K))
  BIC <- rep(.Machine$double.xmax,length(K))
  gMDL <- rep(.Machine$double.xmax,length(K))
  traceS <- rep(.Machine$double.xmax,length(K))
  sigma2list <- rep(.Machine$double.xmax,length(K))
  valpr0 <- eigenvaluesS1[-(1:ddlmini)]
  ImLambda <- rep(0,n)
  ImLambda[-(1:ddlmini)] <- (1-valpr0)
  ImLambdak <- rep(1,n)
  kcourant <- K[1]
  if (kcourant>1) {
    for (j in 1:(kcourant-1)){
      ImLambdak <- ImLambdak*ImLambda
    }
  }
  for (k in K) {
    while (kcourant!=k) {
      ImLambdak <- ImLambdak*ImLambda
      kcourant <- kcourant+1
    }
    ImLambdakp1 <- ImLambdak*ImLambda
    traceSk <- sum(1-ImLambdakp1)
    traceS[iter] <- traceSk
    sigma2 <- (sum((ImLambdakp1^2)*(tUy^2)))/n
    sigma2list[iter] <- sigma2
    logsigma2 <- log(sum((ImLambdakp1^2)*(tUy^2)))-log(n)
    if ((sigma2<1e-10)|(traceSk>((1-1e-10)*n))) break
    GCV[iter] <-  logsigma2 - 2*log(1-traceSk/n) 
    AICc[iter] <-  logsigma2 + 1+ 2*(traceSk+1)/(n-traceSk-2)
    AIC[iter] <-  logsigma2 + 2*(traceSk)/(n)
    BIC[iter] <-  logsigma2 + log(n)*(traceSk)/(n)
    Sbul <- n*sigma2/(n-traceSk)
    gMDL[iter] <- log(Sbul)+traceSk/n*log((sum(y^2)-n*sigma2)/(traceSk*Sbul))
    if (traceSk>ddlmaxi) break
    iter <- iter+1
  }
  
  return(list(aic=AIC,aicc=AICc,gcv=GCV,bic=BIC,gmdl=gMDL,df=traceS,sigma2=sigma2list,kmax=iter))
}
