REML.saery.MA1 <-
function(X, ydi, D, md, sigma2edi, sigma.0 = 1, MAXITER = 20) {
  theta.f <- 0
  sigma1.f <- sigma.0
  sigma2.f <- sigma.0
  thetavect.f  <- c(sigma1.f,sigma2.f,theta.f)
  Bad <- 0
  FLAG <- 0
  p <- ncol(X)
  a <- list(1:md[1])
  mdcum <- cumsum(md)
  for(d in 2:D)
    a[[d]] <- (mdcum[d-1]+1):mdcum[d]
  yd <- Xd <- list()
  for(d in 1:D) {
    yd[[d]] <- ydi[a[[d]]]
    Xd[[d]] <- X[a[[d]],]
  }
  
  for(ITER in 1:MAXITER){
    Vd.inv <- V1d  <- V2d <- V3d <- VinvV1d <- VinvV2d  <- VinvV3d <- Vinvyd <- VinvXd  <- XtVinvV1dVinvX  <- XtVinvV2dVinvX <- XtVinvV3dVinvX <-
      VinvV1dVinvV1d <-  VinvV2dVinvV2d <-  VinvV3dVinvV3d <- XtVinvV1dVinvV1dVinvX <- VinvV1dVinvV2d <- VinvV1dVinvV3d <- VinvV2dVinvV3d <-
      XtVinvV2dVinvV2dVinvX <-  XtVinvV3dVinvV3dVinvX<- XtVinvV1dVinvV2dVinvX <- XtVinvV1dVinvV3dVinvX<- XtVinvV2dVinvV3dVinvX <- XtVinvV3dVinvV3dVinvX <-list()
    Q.inv <- matrix(0, nrow=p, ncol=p)
    tr.VinvV1d  <- tr.VinvV2d <- tr.VinvV3d  <- tr.VinvV1dVinvV1d  <- tr.VinvV2dVinvV2d <- tr.VinvV3dVinvV3d <- tr.VinvV1dVinvV2d  <- tr.VinvV1dVinvV3d <- tr.VinvV2dVinvV3d <- ytVinvX <-             		 ytVinvV1dVinvy  <- SumXtVinvV1dVinvX <-  ytVinvV1dVinvX  <- ytVinvV2dVinvy  <- ytVinvV2dVinvX  <- SumXtVinvV2dVinvX  <- ytVinvV3dVinvy  <- ytVinvV3dVinvX  <- SumXtVinvV3dVinvX  <- 0
    for(d in 1:D) { 
      uno<-rep(1,md[d])
      V1d[[d]] <- uno%*%t(uno)
      Omegad <- matrix(0,nrow=md[d],ncol=md[d])
      for(i in 2:md[d]){
        Omegad[i,i-1]=-theta.f
        Omegad[i-1,i]=-theta.f
      }
      diag(Omegad) <- 1+theta.f^2
      V2d[[d]] <- Omegad
      OmegadPrima <- matrix(0,nrow=md[d],ncol=md[d])
      for(i in 2:md[d]){
        OmegadPrima[i,i-1]=-1
        OmegadPrima[i-1,i]=-1
      }
      diag(OmegadPrima) <- 2*theta.f
      V3d[[d]] <- sigma2.f*OmegadPrima
      
      Vd <- (sigma1.f *(uno%*%t(uno)) + sigma2.f * Omegad + diag(sigma2edi[a[[d]]]))
      
      if(abs(det(Vd))<0.000000001 || abs(det(Vd))>10000000000) {
        FLAG <- 1
        Bad <- Bad+1
        break
      }
      Vd.inv[[d]] <- solve(Vd)                                    
      Vinvyd[[d]] <- Vd.inv[[d]]%*%yd[[d]]                       
      VinvXd[[d]] <- Vd.inv[[d]]%*%Xd[[d]]                       
      Q.inv <- Q.inv + t(Xd[[d]])%*%VinvXd[[d]]                              
      
      VinvV1d[[d]] <- Vd.inv[[d]]%*%V1d[[d]]
      tr.VinvV1d <- tr.VinvV1d + sum(diag(VinvV1d[[d]]))
      XtVinvV1dVinvX[[d]]<- t(VinvXd[[d]])%*%V1d[[d]]%*%VinvXd[[d]]
      ytVinvX <- ytVinvX + t(yd[[d]])%*%VinvXd[[d]]
      ytVinvV1dVinvy <- ytVinvV1dVinvy + t(Vinvyd[[d]])%*%V1d[[d]]%*%Vinvyd[[d]]
      ytVinvV1dVinvX <- ytVinvV1dVinvX + t(Vinvyd[[d]])%*%V1d[[d]]%*%VinvXd[[d]]
      SumXtVinvV1dVinvX <- SumXtVinvV1dVinvX + XtVinvV1dVinvX[[d]]
      
      VinvV2d[[d]] <- Vd.inv[[d]]%*%V2d[[d]]
      tr.VinvV2d <- tr.VinvV2d + sum(diag(VinvV2d[[d]]))
      XtVinvV2dVinvX[[d]] <- t(VinvXd[[d]])%*%V2d[[d]]%*%VinvXd[[d]]
      ytVinvV2dVinvy <- ytVinvV2dVinvy + t(Vinvyd[[d]])%*%V2d[[d]]%*%Vinvyd[[d]]
      ytVinvV2dVinvX <- ytVinvV2dVinvX + t(Vinvyd[[d]])%*%V2d[[d]]%*%VinvXd[[d]]
      SumXtVinvV2dVinvX <- SumXtVinvV2dVinvX + XtVinvV2dVinvX[[d]]
      
      VinvV3d[[d]] <- Vd.inv[[d]]%*%V3d[[d]]
      tr.VinvV3d <- tr.VinvV3d + sum(diag(VinvV3d[[d]]))
      XtVinvV3dVinvX[[d]] <- t(VinvXd[[d]])%*%V3d[[d]]%*%VinvXd[[d]]
      ytVinvV3dVinvy <- ytVinvV3dVinvy + t(Vinvyd[[d]])%*%V3d[[d]]%*%Vinvyd[[d]]
      ytVinvV3dVinvX <- ytVinvV3dVinvX + t(Vinvyd[[d]])%*%V3d[[d]]%*%VinvXd[[d]]
      SumXtVinvV3dVinvX <- SumXtVinvV3dVinvX + XtVinvV3dVinvX[[d]]
      
      VinvV1dVinvV1d[[d]] <- VinvV1d[[d]]%*%VinvV1d[[d]]    
      tr.VinvV1dVinvV1d <- tr.VinvV1dVinvV1d + sum(diag(VinvV1dVinvV1d[[d]]))
      XtVinvV1dVinvV1dVinvX[[d]] <- t(VinvXd[[d]])%*%V1d[[d]]%*%VinvV1d[[d]]%*%VinvXd[[d]]    
      VinvV2dVinvV2d[[d]] <- VinvV2d[[d]]%*%VinvV2d[[d]]     
      tr.VinvV2dVinvV2d <- tr.VinvV2dVinvV2d + sum(diag(VinvV2dVinvV2d[[d]]))
      XtVinvV2dVinvV2dVinvX[[d]] <- t(VinvXd[[d]])%*%V2d[[d]]%*%VinvV2d[[d]]%*%VinvXd[[d]]      
      VinvV3dVinvV3d[[d]] <- VinvV3d[[d]]%*%VinvV3d[[d]]    
      tr.VinvV3dVinvV3d <- tr.VinvV3dVinvV3d + sum(diag(VinvV3dVinvV3d[[d]]))
      XtVinvV3dVinvV3dVinvX[[d]] <- t(VinvXd[[d]])%*%V3d[[d]]%*%VinvV3d[[d]]%*%VinvXd[[d]]     
      VinvV1dVinvV2d[[d]] <- VinvV1d[[d]]%*%VinvV2d[[d]]    
      tr.VinvV1dVinvV2d <- tr.VinvV1dVinvV2d + sum(diag(VinvV1dVinvV2d[[d]]))
      XtVinvV1dVinvV2dVinvX[[d]] <- t(VinvXd[[d]])%*%V1d[[d]]%*%VinvV2d[[d]]%*%VinvXd[[d]]    
      VinvV1dVinvV3d[[d]] <- VinvV1d[[d]]%*%VinvV3d[[d]]      
      tr.VinvV1dVinvV3d <- tr.VinvV1dVinvV3d + sum(diag(VinvV1dVinvV3d[[d]]))
      XtVinvV1dVinvV3dVinvX[[d]] <- t(VinvXd[[d]])%*%V1d[[d]]%*%VinvV3d[[d]]%*%VinvXd[[d]]    
      VinvV2dVinvV3d[[d]] <- VinvV2d[[d]]%*%VinvV3d[[d]]       
      tr.VinvV2dVinvV3d <- tr.VinvV2dVinvV3d + sum(diag(VinvV2dVinvV3d[[d]]))
      XtVinvV2dVinvV3dVinvX[[d]] <- t(VinvXd[[d]])%*%V2d[[d]]%*%VinvV3d[[d]]%*%VinvXd[[d]]
    }
    if(FLAG==1){
      FLAG<-0
      ITER==MAXITER 
      break
    }
    Q <- solve(Q.inv)
    tr.XtVinvV1dVinvXQ <- tr.XtVinvV2dVinvXQ  <- tr.XtVinvV3dVinvXQ <- tr.XtVinvV1dVinvV1dVinvXQ  <- tr.XtVinvV2dVinvV2dVinvXQ <-
      tr.XtVinvV3dVinvV3dVinvXQ <- tr.XtVinvV1dVinvV2dVinvXQ <- tr.XtVinvV1dVinvV3dVinvXQ <- tr.XtVinvV2dVinvV3dVinvXQ <- 
      XtVinvV1dVinvXQ <- XtVinvV2dVinvXQ <- XtVinvV3dVinvXQ <- 0
    for(d in 1:D){
      tr.XtVinvV1dVinvXQ <- tr.XtVinvV1dVinvXQ + sum(diag(XtVinvV1dVinvX[[d]]%*%Q))
      tr.XtVinvV2dVinvXQ <- tr.XtVinvV2dVinvXQ + sum(diag(XtVinvV2dVinvX[[d]]%*%Q))
      tr.XtVinvV3dVinvXQ <- tr.XtVinvV3dVinvXQ + sum(diag(XtVinvV3dVinvX[[d]]%*%Q))
      tr.XtVinvV1dVinvV1dVinvXQ <- tr.XtVinvV1dVinvV1dVinvXQ  + sum(diag(XtVinvV1dVinvV1dVinvX[[d]]%*%Q))
      XtVinvV1dVinvXQ  <- XtVinvV1dVinvXQ + XtVinvV1dVinvX[[d]]%*%Q
      tr.XtVinvV2dVinvV2dVinvXQ <- tr.XtVinvV2dVinvV2dVinvXQ + sum(diag(XtVinvV2dVinvV2dVinvX[[d]]%*%Q))
      XtVinvV2dVinvXQ  <- XtVinvV2dVinvXQ + XtVinvV2dVinvX[[d]]%*%Q      
      tr.XtVinvV3dVinvV3dVinvXQ <- tr.XtVinvV3dVinvV3dVinvXQ + sum(diag(XtVinvV3dVinvV3dVinvX[[d]]%*%Q))
      XtVinvV3dVinvXQ  <- XtVinvV3dVinvXQ + XtVinvV3dVinvX[[d]]%*%Q      
      tr.XtVinvV1dVinvV2dVinvXQ <- tr.XtVinvV1dVinvV2dVinvXQ + sum(diag(XtVinvV1dVinvV2dVinvX[[d]]%*%Q))
      tr.XtVinvV1dVinvV3dVinvXQ <- tr.XtVinvV1dVinvV3dVinvXQ + sum(diag(XtVinvV1dVinvV3dVinvX[[d]]%*%Q))
      tr.XtVinvV2dVinvV3dVinvXQ <- tr.XtVinvV2dVinvV3dVinvXQ + sum(diag(XtVinvV2dVinvV3dVinvX[[d]]%*%Q))
    }
    tr.XtVinvV1dVinvXQXtVinvV1dVinvXQ <- sum(diag(XtVinvV1dVinvXQ%*%XtVinvV1dVinvXQ))
    tr.XtVinvV2dVinvXQXtVinvV2dVinvXQ <- sum(diag(XtVinvV2dVinvXQ%*%XtVinvV2dVinvXQ))
    tr.XtVinvV3dVinvXQXtVinvV3dVinvXQ <- sum(diag(XtVinvV3dVinvXQ%*%XtVinvV3dVinvXQ))
    tr.XtVinvV1dVinvXQXtVinvV2dVinvXQ  <- sum(diag(XtVinvV1dVinvXQ%*%XtVinvV2dVinvXQ))
    tr.XtVinvV1dVinvXQXtVinvV3dVinvXQ  <- sum(diag(XtVinvV1dVinvXQ%*%XtVinvV3dVinvXQ))
    tr.XtVinvV2dVinvXQXtVinvV3dVinvXQ  <- sum(diag(XtVinvV2dVinvXQ%*%XtVinvV3dVinvXQ))
    tr.PV1 <- tr.VinvV1d - tr.XtVinvV1dVinvXQ 
    tr.PV2 <- tr.VinvV2d - tr.XtVinvV2dVinvXQ
    tr.PV3 <- tr.VinvV3d - tr.XtVinvV3dVinvXQ
    tr.PV1PV1  <- tr.VinvV1dVinvV1d - 2*tr.XtVinvV1dVinvV1dVinvXQ + tr.XtVinvV1dVinvXQXtVinvV1dVinvXQ
    tr.PV2PV2  <- tr.VinvV2dVinvV2d - 2*tr.XtVinvV2dVinvV2dVinvXQ + tr.XtVinvV2dVinvXQXtVinvV2dVinvXQ
    tr.PV3PV3  <- tr.VinvV3dVinvV3d - 2*tr.XtVinvV3dVinvV3dVinvXQ + tr.XtVinvV3dVinvXQXtVinvV3dVinvXQ
    tr.PV1PV2  <- tr.VinvV1dVinvV2d - 2*tr.XtVinvV1dVinvV2dVinvXQ + tr.XtVinvV1dVinvXQXtVinvV2dVinvXQ
    tr.PV1PV3  <- tr.VinvV1dVinvV3d - 2*tr.XtVinvV1dVinvV3dVinvXQ + tr.XtVinvV1dVinvXQXtVinvV3dVinvXQ
    tr.PV2PV3  <- tr.VinvV2dVinvV3d - 2*tr.XtVinvV2dVinvV3dVinvXQ + tr.XtVinvV2dVinvXQXtVinvV3dVinvXQ
    ytPV1Py  <- ytVinvV1dVinvy - 2*ytVinvV1dVinvX%*%Q%*%t(ytVinvX)+ytVinvX%*%Q%*%SumXtVinvV1dVinvX%*%Q%*%t(ytVinvX)
    ytPV2Py  <- ytVinvV2dVinvy - 2*ytVinvV2dVinvX%*%Q%*%t(ytVinvX)+ytVinvX%*%Q%*%SumXtVinvV2dVinvX%*%Q%*%t(ytVinvX)
    ytPV3Py  <- ytVinvV3dVinvy - 2*ytVinvV3dVinvX%*%Q%*%t(ytVinvX)+ytVinvX%*%Q%*%SumXtVinvV3dVinvX%*%Q%*%t(ytVinvX)
    S1 <- -0.5*tr.PV1 + 0.5*ytPV1Py
    S2 <- -0.5*tr.PV2 + 0.5*ytPV2Py
    S3 <- -0.5*tr.PV3 + 0.5*ytPV3Py
    F11  <- 0.5*tr.PV1PV1
    F22  <- 0.5*tr.PV2PV2
    F33  <- 0.5*tr.PV3PV3
    F12  <- 0.5*tr.PV1PV2
    F13  <- 0.5*tr.PV1PV3
    F23  <- 0.5*tr.PV2PV3
    Ssig <- c(S1,S2,S3)
    Fsig <- matrix(c(F11,F12,F13,F12,F22,F23,F13,F23,F33),ncol=3)
    if(abs(det(Fsig))< 0.000000001 || sum(abs(Fsig))>10000000000) {
      ITER <- MAXITER
      Bad <- Bad+1
      break
    }
    Fsig.inv <- solve(Fsig)
    dif <- Fsig.inv%*%Ssig
    thetavect.f <- thetavect.f + dif
    theta.f <- thetavect.f[3,1]
    sigma1.f <- thetavect.f[1,1]
    sigma2.f <- thetavect.f[2,1]
    if(identical(as.numeric(abs(dif)<0.00001),rep(1,3))) 
      break
  }
  if(sigma1.f<0 || sigma2.f<0 || theta.f< -1 || theta.f>1) {
    ITER <- MAXITER
    Bad <- Bad+1
  }
  
  return(list(as.vector(thetavect.f), Fsig, ITER, Bad, Q))
  
}
