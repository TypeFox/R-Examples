iterchoiceS1lrcv <- function(X,y,lambda,rank,bs,listvarx,df,ddlmini,ntest,ntrain,Kfold,type,npermut,seed,Kmin,Kmax,criterion,m,s,fraction) {
  choixssecv2 <- function(k,sel,y,lambdalist,rank,S1,valpr,SSx,tUy,ddlmin,index0) {
    sse <- 0
    for (j in 1:length(sel)) {
      if (attr(sel,"type")=="timeseries") {
        nj <- length(y[-(sel[[j]][1]:n)])
      } else {
        nj <- length(y[-sel[[j]]])
    }
      sommerk <- rev(sumvalpr(k,nj,rev(valpr[[j]]),nj-index0[j]+1,nj-ddlmin[j]+1))
      beta <- S1[[j]]$vectors%*%(sommerk[1:rank]*tUy[[j]][1:rank])
      beta <- t((1-valpr[[j]][1:rank])*t(S1[[j]]$Rm1U)) %*%t(S1[[j]]$vectors) %*% beta
      Yrescv <- as.vector(SSx[[j]]%*%beta)
      sse <- sse+sum((y[sel[[j]]]-Yrescv)^2)
    }
    return(sse)
  }
  choixsapcv2 <- function(k,sel,y,lambdalist,rank,S1,valpr,SSx,tUy,ddlmin,index0) {
    sap <- 0
    for (j in 1:length(sel)) {
      if (attr(sel,"type")=="timeseries") {
        nj <- length(y[-(sel[[j]][1]:n)])
      } else {
        nj <- length(y[-sel[[j]]])
      }
      sommerk <- rev(sumvalpr(k,nj,rev(valpr[[j]]),nj-index0[j]+1,nj-ddlmin[j]+1))
      beta <- S1[[j]]$vectors%*%(sommerk[1:rank]*tUy[[j]][1:rank])
      beta <- t((1-valpr[[j]][1:rank])*t(S1[[j]]$Rm1U)) %*%t(S1[[j]]$vectors) %*% beta
      Yrescv <- as.vector(SSx[[j]]%*%beta)
      sap <- sap+sum(abs((y[sel[[j]]]-Yrescv)/y[sel[[j]]]))
    }
    return(sap)
  }
  
  n <- nrow(X)
  sel <- cvobs(n,ntest,ntrain,Kfold,type,npermut,seed)
  res <- list(minimum=0,objective=.Machine$double.xmax)
  objectif <- .Machine$double.xmax
  fraction <- sort(unique(c(fraction[(fraction<Kmax)&(fraction>Kmin)],Kmin,Kmax)))
  dep <- length(fraction)
  SSx <- as.list(rep(0,length(sel)))
  S1 <- valpr <- tUy <-lambdalist <- SSx
  ddlmin <- index0 <-  rep(0,length(sel))
  for (j in 1:length(sel)) {
    if (attr(sel,"type")=="timeseries") {
      XA <- X[-(sel[[j]][1]:n),,drop=FALSE]
      YA <- y[-(sel[[j]][1]:n)]
    } else  {
      XA <- X[-sel[[j]],,drop=FALSE]
      YA <- y[-sel[[j]]]
    }  
    if (is.null(lambda)&(!is.null(df))) {
      lambdalist[[j]] <- lambdachoicelr(XA,ddlmini*df,m,s,rank,itermax=100,bs,listvarx)
    } else lambdalist[[j]] <- lambda
    S1[[j]] <- lrsmoother(XA, bs,listvarx,lambdalist[[j]],m,s,rank) 
    tUy[[j]] <- as.vector(crossprod(S1[[j]]$vectors, YA))
    if (any(zapsmall(S1[[j]]$values-1,digits=9)==0)) {
      ddlmin[j] <-  sum(zapsmall(S1[[j]]$values-1,digits=9)==0)
    } else ddlmin[j] <- NA
    if (any(zapsmall(S1[[j]]$values,digits=9)==0)) {
      index0[j] <-  which(zapsmall(S1[[j]]$values,digits=9)==0)[1]
    } else index0[j] <- rank+1
    valpr[[j]] <- 1-c(S1[[j]]$values,rep(0,nrow(XA)-rank))
    valpr[[j]][valpr[[j]]<0] <- 0
    SSx[[j]] <- PredictMat(S1[[j]]$smoothobject,data.frame(X[sel[[j]],,drop=FALSE]))
    S1[[j]]$smoothobject <- NULL
  }
  if (criterion=="rmse") {
    repeat {
      mini <- fraction[dep-1]
      objectif <- choixssecv2(mini,sel=sel,y=y,lambdalist=lambdalist,rank=rank,S1=S1,valpr=valpr,SSx=SSx,tUy=tUy,ddlmin=ddlmin,index0=index0)
      if (objectif<.Machine$double.xmax) {
        bb <- fraction[dep]
        break 
      }
      dep <- dep-1
      if (dep==1) stop(paste("decrease Kmax below",fraction[2]))
    }
    phi <- (sqrt(5) - 1)/2
    repeat {
      x1 <- mini + (1-phi)*(bb-mini)
      objectif <- choixssecv2(x1,sel=sel,y=y,lambdalist=lambdalist,rank=rank,S1=S1,valpr=valpr,SSx=SSx,tUy=tUy,ddlmin=ddlmin,index0=index0)
      if (objectif<.Machine$double.xmax) break else {
        bb <- x1
      }
    }
    fraction <- c(fraction[1:(dep-1)],bb)
    for (i in 1:(length(fraction)-1)) {
      res1 <- optimize(choixssecv2,lower=fraction[i],upper=fraction[i+1],tol=0.5,sel=sel,y=y,lambdalist=lambdalist,rank=rank,S1=S1,valpr=valpr,SSx=SSx,tUy=tUy,ddlmin=ddlmin,index0=index0)
      if (res1$objective<res$objective) res <- res1
    }
  } else {
    repeat {
      mini <- fraction[dep-1]
      objectif <- choixsapcv2(mini,sel=sel,y=y,lambdalist=lambdalist,rank=rank,S1=S1,valpr=valpr,SSx=SSx,tUy=tUy,ddlmin=ddlmin,index0=index0)
      if (objectif<.Machine$double.xmax) {
        bb <- fraction[dep]
        break 
      }
      dep <- dep-1
      if (dep==1) stop(paste("decrease Kmax below",fraction[2]))
    }
    phi <- (sqrt(5) - 1)/2
    repeat {
      x1 <- mini + (1-phi)*(bb-mini)
      objectif <- choixsapcv2(x1,sel=sel,y=y,lambdalist=lambdalist,rank=rank,S1=S1,valpr=valpr,SSx=SSx,tUy=tUy,ddlmin=ddlmin,index0=index0)
      if (objectif<.Machine$double.xmax) break else {
        bb <- x1
      }
    }
    fraction <- c(fraction[1:(dep-1)],bb)
    for (i in 1:(length(fraction)-1)) {
      res1 <- optimize(choixsapcv2,lower=fraction[i],upper=fraction[i+1],tol=0.5,sel=sel,y=y,lambdalist=lambdalist,rank=rank,S1=S1,valpr=valpr,SSx=SSx,tUy=tUy,ddlmin=ddlmin,index0=index0)
      if (res1$objective<res$objective) res <- res1
    }
  }
  return(list(iter=round(res$minimum),objective=res$objective/sum(unlist(lapply(sel,length)))))
}


