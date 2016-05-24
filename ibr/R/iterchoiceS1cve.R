iterchoiceS1cve <- function(X,y,lambda,df,ddlmini,ntest,ntrain,Kfold,type,npermut,seed,Kmin,Kmax,m,s) {
  iter <- 1
  MSE <- MAP <- rep(.Machine$double.xmax,length(Kmax-Kmin+1))
  n <- nrow(X)
  sel <- cvobs(n,ntest,ntrain,Kfold,type,npermut,seed)
  U <- as.list(rep(0,length(sel)))
  sumvalprk <- lambdalist <- S1 <- valpr <- tUy <- Sp <- preprod <- SSx <- U
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
      lambdalist[[j]] <- lambdachoice(XA,ddlmini*df,m=m,s=s,itermax=100)
    } else lambdalist[[j]] <- lambda
    nj[j] <- length(YA)
    S1[[j]] <- dssmoother(XA, YA,lambda=lambdalist[[j]],m=m,s=s) 
    vp1.S1 <- eigen(S1[[j]]$H,symmetric=TRUE)
    U[[j]] <- vp1.S1$vect
    tUy[[j]] <- as.vector(crossprod(U[[j]], YA))
    valpr0 <- vp1.S1$values[-(1:ddlmini)]
    valpr[[j]] <- rep(0,nj[j])
    valpr[[j]][-(1:ddlmini)] <- (1-valpr0)
    sumvalprk[[j]] <- rep(1,nj[j])
    qrSgu <- qr(S1[[j]]$Sgu)
    F2 <- qr.Q(qrSgu,complete=TRUE)[,-(1:ncol(S1[[j]]$Sgu))]
    ainv <- t(F2)%*%S1[[j]]$Qgu%*%F2
    diag(ainv) <- diag(ainv)+lambdalist[[j]]
    Sp[[j]] <- -lambda*F2%*%(solve(ainv))%*%t(F2)
    SSx[[j]] <- dsSx(X=XA,X[sel[[j]],,drop=FALSE],m,s)
    preprod[[j]] <- solve(qr.R(qrSgu))%*%(t(qr.Q(qrSgu)))
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
      beta <- U[[j]]%*%(sumvalprk[[j]]*tUy[[j]])
      cgubeta <- as.vector(Sp[[j]]%*%beta)/(-lambdalist[[j]])
      dgubeta <- preprod[[j]]%*%(as.matrix(beta)-(S1[[j]])$Qgu%*%cgubeta)
      Yrescv <- as.vector((SSx[[j]])$Sgu%*%dgubeta+(SSx[[j]])$Qgu%*%cgubeta)
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

