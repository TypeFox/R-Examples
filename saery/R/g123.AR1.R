g123.AR1 <-
function(X, D, md, sigma2edi, sigmau1, sigmau2, rho, Fsig) {
  
  p <- ncol(X)
  a <- list(1:md[1])
  mdcum <- cumsum(md)
  for(d in 2:D)
    a[[d]] <- (mdcum[d-1]+1):mdcum[d]
  
  Xd <- Vd.inv <- Ved.inv <- Omegad  <- OmegadPrima <- T11d <- T12d <- T22d <- VinvOmega <- OmegaVinvOmega <-VedinvXd <- auxq <- unotunoVinv <-
    VinvOmegaPrima  <-  OmegaVinvOmegaPrima <-  OmegaPrimaVinvOmegaPrima <-g1.a  <- g2.a  <- q11  <- q12 <- q13 <- q22 <- q23 <- q33 <- list()
  Q.inv <- matrix(0, nrow=p, ncol=p)
  
  for(d in 1:D) {
    uno<-rep(1,md[d])
    unotuno<-uno%*%t(uno)
    
    Omegad[[d]]<-matrix(0,nrow=md[d],ncol=md[d])
    Omegad[[d]][lower.tri(Omegad[[d]])]<-rho^sequence((md[d]-1):1)
    Omegad[[d]]<-Omegad[[d]]+t(Omegad[[d]])
    diag(Omegad[[d]])<-1
    Omegad[[d]] <- (1/(1-rho^2))*Omegad[[d]]
    
    OmegadPrima[[d]]<-matrix(0,nrow=md[d],ncol=md[d])
    OmegadPrima[[d]][lower.tri(OmegadPrima[[d]])]<-sequence((md[d]-1):1)*rho^(sequence((md[d]-1):1)-1)
    OmegadPrima[[d]]<-OmegadPrima[[d]]+t(OmegadPrima[[d]])
    OmegadPrima[[d]]<- (1/(1-rho^2))*OmegadPrima[[d]]
    OmegadPrima[[d]] <- OmegadPrima[[d]] + (2*rho/(1-rho^2))*Omegad[[d]]
    
    Xd[[d]] <- X[a[[d]],]
    Ved <- diag(sigma2edi[a[[d]]])
    Vd <- (sigmau1 *(uno%*%t(uno)) + sigmau2 * Omegad[[d]] +  Ved)
    Vd.inv[[d]] <- solve(Vd)
    Ved.inv[[d]] <- solve(Ved)
    VinvOmega[[d]] <-  Vd.inv[[d]]%*%  Omegad[[d]]
    OmegaVinvOmega[[d]] <- t(VinvOmega[[d]])%*%Omegad[[d]]
    VedinvXd[[d]] <- Ved.inv[[d]]%*%Xd[[d]]
    VinvOmegaPrima[[d]] <-  Vd.inv[[d]]%*%OmegadPrima[[d]]
    OmegaVinvOmegaPrima[[d]] <- Omegad[[d]]%*%Vd.inv[[d]]%*%OmegadPrima[[d]]
    OmegaPrimaVinvOmegaPrima[[d]] <- OmegadPrima[[d]]%*%Vd.inv[[d]]%*%OmegadPrima[[d]]
    
    T11d[[d]] <- sigmau1-sigmau1^2*t(uno)%*%Vd.inv[[d]]%*%uno
    T12d[[d]] <- -sigmau1*sigmau2*t(uno)%*%VinvOmega[[d]]
    T22d[[d]] <- sigmau2*Omegad[[d]] -sigmau2^2*OmegaVinvOmega[[d]] 
    
    g1.a[[d]] <- uno%*%T11d[[d]]%*%t(uno) + uno%*%T12d[[d]] + t(T12d[[d]])%*%t(uno) + T22d[[d]]
    g2.a[[d]] <- Xd[[d]] - (uno%*%T11d[[d]]%*%t(uno) + uno%*%T12d[[d]] + t(T12d[[d]])%*%t(uno) + T22d[[d]])%*% VedinvXd[[d]]
    
    auxq[[d]] <- sigmau1*unotuno + sigmau2* Omegad[[d]]
    unotunoVinv[[d]] <- unotuno%*%Vd.inv[[d]]
    q11[[d]] <- unotunoVinv[[d]]%*%unotuno - 2*unotunoVinv[[d]]%*%unotunoVinv[[d]]%*%auxq[[d]]+auxq[[d]]%*%Vd.inv[[d]]%*%unotunoVinv[[d]]%*%unotunoVinv[[d]]%*%auxq[[d]]
    q12[[d]]  <- unotuno%*%VinvOmega[[d]]-unotunoVinv[[d]]%*%t(VinvOmega[[d]])%*%auxq[[d]] - auxq[[d]]%*%t(unotunoVinv[[d]])%*%VinvOmega[[d]] + auxq[[d]]%*%t(unotunoVinv[[d]])%*%VinvOmega[[d]]%*%Vd.inv[[d]]%*%auxq[[d]]
    q13[[d]] <- sigmau2*unotuno%*%VinvOmegaPrima[[d]] - sigmau2*unotunoVinv[[d]]%*%t(VinvOmegaPrima[[d]])%*%auxq[[d]] - sigmau2*auxq[[d]]%*%t(unotunoVinv[[d]])%*%VinvOmegaPrima[[d]] + sigmau2*auxq[[d]]%*%Vd.inv[[d]]%*%unotunoVinv[[d]]%*%t(VinvOmegaPrima[[d]]%*%auxq[[d]])
    q22[[d]]  <- OmegaVinvOmega[[d]] - 2*OmegaVinvOmega[[d]]%*%Vd.inv[[d]]%*%auxq[[d]] + auxq[[d]]%*%Vd.inv[[d]]%*%OmegaVinvOmega[[d]]%*%Vd.inv[[d]]%*%auxq[[d]]    
    q23[[d]] <- sigmau2*OmegaVinvOmegaPrima[[d]]  - sigmau2*OmegaVinvOmegaPrima[[d]] %*%Vd.inv[[d]]%*%auxq[[d]] - sigmau2*auxq[[d]]%*%VinvOmega[[d]]%*%VinvOmegaPrima[[d]] + sigmau2*auxq[[d]]%*%Vd.inv[[d]]%*%OmegaVinvOmegaPrima[[d]]%*%Vd.inv[[d]]%*%auxq[[d]]	
    q33[[d]] <- sigmau2^2*OmegaPrimaVinvOmegaPrima[[d]] - 2*sigmau2^2*OmegaPrimaVinvOmegaPrima[[d]]%*%Vd.inv[[d]]%*%auxq[[d]] + sigmau2^2*auxq[[d]]%*%Vd.inv[[d]]%*%OmegaPrimaVinvOmegaPrima[[d]]%*%Vd.inv[[d]]%*%auxq[[d]]
    
    Q.inv <- Q.inv + t(Xd[[d]])%*%Vd.inv[[d]]%*%Xd[[d]]
    
  }
  Q <- solve(Q.inv)
  
  g1 <- g2 <- g3 <- list()
  for(d in 1:D){
    g1[[d]] <- diag(g1.a[[d]])
    g2[[d]] <- diag(g2.a[[d]]%*%Q%*%t(g2.a[[d]]))
    q11[[d]] <- diag(q11[[d]])
    q12[[d]] <- diag(q12[[d]])
    q13[[d]] <- diag(q13[[d]])
    q22[[d]] <- diag(q22[[d]])
    q23[[d]] <- diag(q23[[d]])
    q33[[d]] <- diag(q33[[d]])
  }
  for(d in 1:D){
    g3[[d]] <- vector()
    for(i in 1:md[d]){
      g3[[d]][i] <- sum(diag(matrix(c(q11[[d]][i],q12[[d]][i],q13[[d]][i],q12[[d]][i],q22[[d]][i],q23[[d]][i],q13[[d]][i],q23[[d]][i],q33[[d]][i]),nrow=3)%*%solve(Fsig)))
    }
  }
  g1 <- unlist(g1)
  g2 <- unlist(g2)
  g3 <- unlist(g3)
  
  return(list(g1,g2,g3))
  
}
