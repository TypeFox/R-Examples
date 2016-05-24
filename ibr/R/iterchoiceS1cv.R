iterchoiceS1cv <- function(X,y,lambda,df,ddlmini,ntest,ntrain,Kfold,type,npermut,seed,Kmin,Kmax,criterion,m,s,fraction) {
  choixssecv2 <- function(k,sel,y,lambdalist,U,S1,valpr,SSx,tUy,Sp,ddlmin,index0) {
    sse <- 0
    for (j in 1:length(sel)) {
      if (attr(sel,"type")=="timeseries") {
        nj <- length(y[-(sel[[j]][1]:n)])
      } else {
        nj <- length(y[-sel[[j]]])
      }
      sommerk <- rev(sumvalpr(k,nj,rev(valpr[[j]]),nj-index0[j]+1,nj-ddlmin[j]+1))
      beta <- U[[j]]%*%(sommerk*tUy[[j]])
      cgubeta <- as.vector(Sp[[j]]%*%beta)/(-lambdalist[[j]])
      dgubeta <- preprod[[j]]%*%(as.matrix(beta)-(S1[[j]])$Qgu%*%cgubeta)
      Yrescv <- as.vector((SSx[[j]])$Sgu%*%dgubeta+(SSx[[j]])$Qgu%*%cgubeta)
      sse <- sse+sum((y[sel[[j]]]-Yrescv)^2)
    }
    return(sse)
  }
  choixsapcv2 <- function(k,sel,y,lambdalist,U,S1,valpr,SSx,tUy,Sp,ddlmin,index0) {
    sap <- 0
    for (j in 1:length(sel)) {
      if (attr(sel,"type")=="timeseries") {
        nj <- length(y[-(sel[[j]][1]:n)])
      } else {
        nj <- length(y[-sel[[j]]])
      }
      sommerk <- rev(sumvalpr(k,nj,rev(valpr[[j]]),nj-index0[j]+1,nj-ddlmin[j]+1))
      beta <- U[[j]]%*%(sommerk*tUy[[j]])
      cgubeta <- as.vector(Sp[[j]]%*%beta)/(-lambdalist[[j]])
      dgubeta <- preprod[[j]]%*%(as.matrix(beta)-(S1[[j]])$Qgu%*%cgubeta)
      Yrescv <- as.vector((SSx[[j]])$Sgu%*%dgubeta+(SSx[[j]])$Qgu%*%cgubeta)
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
  U <- as.list(rep(0,length(sel)))
  S1 <- valpr <- tUy <- Sp <- lambdalist <- preprod <- SSx <- U
  ddlmin <- index0 <- rep(0,length(sel))
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
    S1[[j]] <- dssmoother(XA, YA,lambda=lambdalist[[j]],m=m,s=s) 
    vp1.S1 <- eigen(S1[[j]]$H,symmetric=TRUE)
    U[[j]] <- vp1.S1$vect
    tUy[[j]] <- as.vector(crossprod(U[[j]], YA))
    if (any(zapsmall(vp1.S1$values-1,digits=9)==0)) {
      ddlmin[j] <-  sum(zapsmall(vp1.S1$values-1,digits=9)==0)
    } else ddlmin[j] <- NA
    if (any(zapsmall(vp1.S1$values,digits=9)==0)) {
      index0[j] <-  which(zapsmall(vp1.S1$values,digits=9)==0)[1]
    } else index0[j] <- NA
    valpr[[j]] <- 1-vp1.S1$values
    qrSgu <- qr(S1[[j]]$Sgu)
    F2 <- qr.Q(qrSgu,complete=TRUE)[,-(1:ncol(S1[[j]]$Sgu))]
    ainv <- t(F2)%*%S1[[j]]$Qgu%*%F2
    diag(ainv) <- diag(ainv)+lambdalist[[j]]
    Sp[[j]] <- -lambda*F2%*%(solve(ainv))%*%t(F2)
    SSx[[j]] <- dsSx(X=XA,X[sel[[j]],,drop=FALSE],m,s)
    preprod[[j]] <- solve(qr.R(qrSgu))%*%(t(qr.Q(qrSgu)))
  }
  if (criterion=="rmse") {
    repeat {
      mini <- fraction[dep-1]
      objectif <- choixssecv2(mini,sel=sel,y=y,lambdalist=lambdalist,U=U,S1=S1,valpr=valpr,SSx=SSx,tUy=tUy,Sp=Sp,ddlmin=ddlmin,index0=index0)
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
      objectif <- choixssecv2(x1,sel=sel,y=y,lambdalist=lambdalist,U=U,S1=S1,valpr=valpr,SSx=SSx,tUy=tUy,Sp=Sp,ddlmin=ddlmin,index0=index0)
      if (objectif<.Machine$double.xmax) break else {
        bb <- x1
      }
    }
    fraction <- c(fraction[1:(dep-1)],bb)
    for (i in 1:(length(fraction)-1)) {
      res1 <- optimize(choixssecv2,lower=fraction[i],upper=fraction[i+1],tol=0.5,sel=sel,y=y,lambdalist=lambdalist,U=U,S1=S1,valpr=valpr,SSx=SSx,tUy=tUy,Sp=Sp,ddlmin=ddlmin,index0=index0)
      if (res1$objective<res$objective) res <- res1
    }
  } else {
    repeat {
      mini <- fraction[dep-1]
      objectif <- choixsapcv2(mini,sel=sel,y=y,lambdalist=lambdalist,U=U,S1=S1,valpr=valpr,SSx=SSx,tUy=tUy,Sp=Sp,ddlmin=ddlmin,index0=index0)
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
      objectif <- choixsapcv2(x1,sel=sel,y=y,lambdalist=lambdalist,U=U,S1=S1,valpr=valpr,SSx=SSx,tUy=tUy,Sp=Sp,ddlmin=ddlmin,index0=index0)
      if (objectif<.Machine$double.xmax) break else {
        bb <- x1
      }
    }
    fraction <- c(fraction[1:(dep-1)],bb)
    for (i in 1:(length(fraction)-1)) {
      res1 <- optimize(choixsapcv2,lower=fraction[i],upper=fraction[i+1],tol=0.5,sel=sel,y=y,lambdalist=lambdalist,U=U,S1=S1,valpr=valpr,SSx=SSx,tUy=tUy,Sp=Sp,ddlmin=ddlmin,index0=index0)
      if (res1$objective<res$objective) res <- res1
    }
  }
  return(list(iter=round(res$minimum),objective=res$objective/sum(unlist(lapply(sel,length)))))
}


