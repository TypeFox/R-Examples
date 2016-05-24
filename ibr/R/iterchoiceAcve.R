iterchoiceAcve <- function(X,y,bx,df=NULL,kernelx,ddlmini,ntest,ntrain,Kfold,type,npermut,seed,Kmin,Kmax) {
  iter <- 1
  MSE <- MAP <- rep(.Machine$double.xmax,length(Kmax-Kmin+1))
  n <- nrow(X)
  sel <- cvobs(n,ntest,ntrain,Kfold,type,npermut,seed)
    tPADmdemiY <- as.list(rep(0,length(sel)))
    sumvalprk <- valpr <- DdemiPA <- SSx <- tPADmdemiY
    nj <- rep(0,length(sel))
     for (j in 1:length(sel)) {
      if (attr(sel,"type")=="timeseries") {
        XA <- X[-(sel[[j]][1]:n),,drop=FALSE]
        YA <- y[-(sel[[j]][1]:n)]
      } else  {
        XA <- X[-sel[[j]],,drop=FALSE]
        YA <- y[-sel[[j]]]
      }  
       nj[j] <- length(YA)
      if (is.null(bx)&(!is.null(df))) {
        bx <- bwchoice(X=XA,objectif=df,kernelx,itermax=100)
      }
      listeA <- calcA(X=XA,bx=bx,kernelx=kernelx)
      listeA.eig <- eigen(listeA$A,symmetric=TRUE)
      tPADmdemiY[[j]] <- t(listeA.eig$vectors*(1/listeA$Ddemi))%*%YA
      ddlmini <- sum(zapsmall(listeA.eig$values-1)==0)
      valpr[[j]] <- rep(0,nj[j])
      sumvalprk[[j]] <- rep(1,nj[j])
      if (ddlmini>=1) {
        valpr[[j]][-(1:ddlmini)] <- 1-listeA.eig$values[-(1:ddlmini)]
      } else {
        valpr[[j]] <- 1-listeA.eig$values
      }
      DdemiPA[[j]] <- (listeA$Ddemi*listeA.eig$vectors)
      SSx[[j]] <- kernelSx(kernelx=kernelx,X=XA,bx=bx,X[sel[[j]],,drop=FALSE])
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
        prov1 <- matrix(sumvalprk[[j]]*as.vector(tPADmdemiY[[j]]),nj[j],1)
        Yrescv <- SSx[[j]]%*%DdemiPA[[j]]%*%prov1
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
