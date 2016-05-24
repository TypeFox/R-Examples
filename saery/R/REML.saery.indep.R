REML.saery.indep <-
function(X, ydi, D, md, sigma2edi, sigma1.0 = 1, sigma2.0 = 1, MAXITER = 20) {
  sigma1.f <- sigma1.0
  sigma2.f <- sigma2.0
  theta.f  <- c(sigma1.f,sigma2.f)
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
    Vd.inv <- Vinvyd <- VinvXd <- list()
    VinvV1d <- XtVinvV1dVinvX <- VinvV1dVinvV1d <- XtVinvV1dVinvV1dVinvX <- list()
    tr.VinvV1d <- ytVinvX <- ytVinvV1dVinvy <- ytVinvV1dVinvX <- SumXtVinvV1dVinvX <- tr.VinvV1dVinvV1d <- 0
    VinvV1dVinvV2d <- XtVinvV1dVinvV2dVinvX <- list()
    tr.VinvV1dVinvV2d <- 0
    VinvV2d <- XtVinvV2dVinvX <- VinvV2dVinvV2d <- XtVinvV2dVinvV2dVinvX <- list()
    tr.VinvV2d <- ytVinvV2dVinvy <- ytVinvV2dVinvX <- SumXtVinvV2dVinvX <- tr.VinvV2dVinvV2d <- 0
    Q.inv <- matrix(0, nrow=p, ncol=p)
    for(d in 1:D) {
      uno <- rep(1,md[d])
      V1d <- uno%*%t(uno)
      V2d <- diag(uno)
      Vd <- (sigma1.f*V1d + sigma2.f*V2d + diag(sigma2edi[a[[d]]]))
      if(abs(det(Vd))<0.000000001 || abs(det(Vd))>100000000000) {
        FLAG <- 1
        Bad <- Bad+1
        break
      }
      Vd.inv[[d]] <- solve(Vd)
      Vinvyd[[d]] <- Vd.inv[[d]]%*%yd[[d]]
      VinvXd[[d]] <- Vd.inv[[d]]%*%Xd[[d]]
      Q.inv <- Q.inv + t(Xd[[d]])%*%VinvXd[[d]]
      VinvV1d[[d]] <- Vd.inv[[d]]%*%V1d
      tr.VinvV1d <- tr.VinvV1d + sum(diag(VinvV1d[[d]]))
      XtVinvV1dVinvX[[d]] <- t(VinvXd[[d]])%*%V1d%*%VinvXd[[d]]
      ytVinvX <- ytVinvX + t(yd[[d]])%*%VinvXd[[d]]
      ytVinvV1dVinvy <- ytVinvV1dVinvy + t(Vinvyd[[d]])%*%V1d%*%Vinvyd[[d]]
      ytVinvV1dVinvX <- ytVinvV1dVinvX + t(Vinvyd[[d]])%*%V1d%*%VinvXd[[d]]
      SumXtVinvV1dVinvX <- SumXtVinvV1dVinvX + XtVinvV1dVinvX[[d]]
      VinvV2d[[d]] <- Vd.inv[[d]]%*%V2d
      tr.VinvV2d <-  tr.VinvV2d + sum(diag(VinvV2d[[d]]))
      XtVinvV2dVinvX[[d]] <- t(VinvXd[[d]])%*%V2d%*%VinvXd[[d]]
      ytVinvV2dVinvy <- ytVinvV2dVinvy + t(Vinvyd[[d]])%*%V2d%*%Vinvyd[[d]]
      ytVinvV2dVinvX <- ytVinvV2dVinvX + t(Vinvyd[[d]])%*%V2d%*%VinvXd[[d]]
      SumXtVinvV2dVinvX <- SumXtVinvV2dVinvX + XtVinvV2dVinvX[[d]]
      VinvV1dVinvV1d[[d]] <- VinvV1d[[d]]%*%VinvV1d[[d]]
      tr.VinvV1dVinvV1d <- tr.VinvV1dVinvV1d + sum(diag(VinvV1dVinvV1d[[d]]))
      XtVinvV1dVinvV1dVinvX[[d]]  <- t(Xd[[d]])%*%VinvV1dVinvV1d[[d]]%*%VinvXd[[d]]
      VinvV1dVinvV2d[[d]] <- VinvV1d[[d]]%*%VinvV2d[[d]]
      tr.VinvV1dVinvV2d <- tr.VinvV1dVinvV2d + sum(diag(VinvV1dVinvV2d[[d]]))
      XtVinvV1dVinvV2dVinvX[[d]] <- t(Xd[[d]])%*%VinvV1dVinvV2d[[d]]%*%VinvXd[[d]]
      VinvV2dVinvV2d[[d]] <- VinvV2d[[d]]%*%VinvV2d[[d]]
      tr.VinvV2dVinvV2d <- tr.VinvV2dVinvV2d + sum(diag(VinvV2dVinvV2d[[d]]))
      XtVinvV2dVinvV2dVinvX[[d]] <- t(Xd[[d]])%*%VinvV2dVinvV2d[[d]]%*%VinvXd[[d]]
      if(FLAG==1){
        FLAG <- 0
        ITER <- MAXITER 
        break
      }
    }
    Q <- solve(Q.inv)
    tr.XtVinvV1dVinvXQ <- tr.XtVinvV2dVinvXQ <- 0
    tr.XtVinvV1dVinvV1dVinvXQ <- tr.XtVinvV1dVinvV2dVinvXQ <- tr.XtVinvV2dVinvV2dVinvXQ <- 0
    for(d in 1:D){
      tr.XtVinvV1dVinvXQ <- tr.XtVinvV1dVinvXQ + sum(diag(XtVinvV1dVinvX[[d]]%*%Q))
      tr.XtVinvV2dVinvXQ <- tr.XtVinvV2dVinvXQ + sum(diag(XtVinvV2dVinvX[[d]]%*%Q))
      tr.XtVinvV1dVinvV1dVinvXQ <- tr.XtVinvV1dVinvV1dVinvXQ + sum(diag(XtVinvV1dVinvV1dVinvX[[d]]%*%Q))
      tr.XtVinvV1dVinvV2dVinvXQ <- tr.XtVinvV1dVinvV2dVinvXQ + sum(diag(XtVinvV1dVinvV2dVinvX[[d]]%*%Q))
      tr.XtVinvV2dVinvV2dVinvXQ <- tr.XtVinvV2dVinvV2dVinvXQ + sum(diag(XtVinvV2dVinvV2dVinvX[[d]]%*%Q))
    }
    tr.PV1 <- tr.VinvV1d - tr.XtVinvV1dVinvXQ 
    tr.PV2 <- tr.VinvV2d - tr.XtVinvV2dVinvXQ
    SumXtVinvV1dVinvXQ <- SumXtVinvV1dVinvX%*%Q
    SumXtVinvV2dVinvXQ <- SumXtVinvV2dVinvX%*%Q
    tr.PV1PV1  <- tr.VinvV1dVinvV1d - 2*tr.XtVinvV1dVinvV1dVinvXQ + sum(diag(SumXtVinvV1dVinvXQ%*%SumXtVinvV1dVinvXQ))
    tr.PV1PV2  <- tr.VinvV1dVinvV2d - 2*tr.XtVinvV1dVinvV2dVinvXQ + sum(diag(SumXtVinvV1dVinvXQ%*%SumXtVinvV2dVinvXQ))
    tr.PV2PV2  <- tr.VinvV2dVinvV2d - 2*tr.XtVinvV2dVinvV2dVinvXQ + sum(diag(SumXtVinvV2dVinvXQ%*%SumXtVinvV2dVinvXQ))
    ytVinvXQ <- ytVinvX%*%Q
    ytPV1Py  <- ytVinvV1dVinvy - 2*ytVinvV1dVinvX%*%t(ytVinvXQ) + ytVinvXQ%*%SumXtVinvV1dVinvX%*%t(ytVinvXQ)
    ytPV2Py  <- ytVinvV2dVinvy - 2*ytVinvV2dVinvX%*%t(ytVinvXQ) + ytVinvXQ%*%SumXtVinvV2dVinvX%*%t(ytVinvXQ)
    S1 <- -0.5*tr.PV1 + 0.5*ytPV1Py
    S2 <- -0.5*tr.PV2 + 0.5*ytPV2Py
    F11  <- 0.5*tr.PV1PV1
    F12  <- 0.5*tr.PV1PV2
    F22  <- 0.5*tr.PV2PV2
    Ssig <- c(S1,S2)
    Fsig <- matrix(c(F11,F12,F12,F22),ncol=2)
    Fsig.inv <- solve(Fsig)
    dif <- Fsig.inv%*%Ssig
    theta.f <- theta.f + dif
    sigma1.f <- theta.f[1,1]
    sigma2.f <- theta.f[2,1]      
    if(identical(as.numeric(abs(dif)<0.000001),rep(1,2))){
      break
    }
  }
  if(sigma1.f<0 || sigma2.f<0) {
    ITER <- MAXITER
    Bad <- Bad+1
  }
  
  return(list(as.vector(theta.f), Fsig, ITER, Bad, Q))
  
}
