iterchoiceS1lrcve <- function(X,y,lambda,rank,bs,listvarx,df,ddlmini,ntest,ntrain,Kfold,type,npermut,seed,Kmin,Kmax,m,s) {
  iter <- 1
  MSE <- MAP <- rep(.Machine$double.xmax,length(Kmax-Kmin+1))
  n <- nrow(X)
 sel <- cvobs(n,ntest,ntrain,Kfold,type,npermut,seed)
  SSx <- as.list(rep(0,length(sel)))
  sumvalprk <- lambdalist <- S1 <- valpr <- tUy <- SSx
  nj <- rep(0,length(sel))
  for (j in 1:length(sel)) {
    if (attr(sel,"type")=="timeseries") {
      XA <- X[-(sel[[j]][1]:n),,drop=FALSE]
      YA <- y[-(sel[[j]][1]:n)]
    } else  {
      XA <- X[-sel[[j]],,drop=FALSE]
      YA <- y[-sel[[j]]]
    }  
    if (is.null(lambda)&(!is.null(df))) {
      lambdalist[[j]] <- lambdachoicelr(XA,ddlmini*df,m,s,rank,100,bs,listvarx)
    } else lambdalist[[j]] <- lambda
    nj[j] <- length(YA)
    S1[[j]] <- lrsmoother(XA, bs,listvarx,lambdalist[[j]],m,s,rank)
    tUy[[j]] <- as.vector(crossprod(S1[[j]]$vectors, YA))
    valpr0 <- S1[[j]]$values[-(1:ddlmini)]
    valpr[[j]] <- rep(0,rank)
    valpr[[j]][-(1:ddlmini)] <- (1-valpr0)
    sumvalprk[[j]] <- rep(1,rank)
    SSx[[j]] <- PredictMat(S1[[j]]$smoothobject,data.frame(X[sel[[j]],,drop=FALSE]))
  }
  if (Kmin>1) {
    for (k in 1:(Kmin-1)){
      for (j in 1:length(sel)) {
        sumvalprk[[j]] <- sumvalprk[[j]]+valpr[[j]]^k
      }
    }
  }
  for (k in Kmin:Kmax) {
    sse <- 0
    sap <- 0
    for (j in 1:length(sel)) {
  beta <- S1[[j]]$vectors%*%(sumvalprk[[j]]*tUy[[j]][1:rank])
  beta <- t((1-valpr[[j]])*t(S1[[j]]$Rm1U)) %*%t(S1[[j]]$vectors) %*% beta
       Yrescv <- as.vector(SSx[[j]]%*%beta)
      sse <- sse+sum((y[sel[[j]]]-Yrescv)^2)
      sap <- sap+sum(abs((y[sel[[j]]]-Yrescv)/y[sel[[j]]]))
      sumvalprk[[j]] <- sumvalprk[[j]]+valpr[[j]]^k
    }
    MSE[iter] <- sse/sum(n-nj)
    MAP[iter] <- sap/sum(n-nj)
    iter <- iter+1
  }
  return(list(rmse=sqrt(MSE),map=MAP))
}     

