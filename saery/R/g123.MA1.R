g123.MA1 <-
function(X, D, md, sigma2edi, sigmau1, sigmau2, theta, Fsig) {
  
  p <- ncol(X)
  a <- list(1:md[1])
  mdcum <- cumsum(md)
  for(d in 2:D)
    a[[d]] <- (mdcum[d-1]+1):mdcum[d]
  
  
  Xd <- Vd.inv <- Ved.inv <- Omegad  <- OmegadPrima <- T11d <- T12d <- T22d <- VinvOmega <- OmegaVinvOmega <-VedinvXd <- auxq <- unotunoVinv <-
    VinvOmegaPrima  <-  OmegaVinvOmegaPrima <-  OmegaPrimaVinvOmegaPrima <-g1.a  <- g2.a  <- q11  <- q12 <- q13 <- q22 <- q23 <- q33 <- list()                                            #OmegaVinvOmegadSinvXd 
  Q.inv <- matrix(0, nrow=p, ncol=p)
  
  for(d in 1:D) {
    
    uno<-rep(1,md[d])
    unotuno<-uno%*%t(uno)
    
    #Matriz Omegad y su derivada  
    
    Omegad[[d]] <- matrix(0,nrow=md[d],ncol=md[d])
    for(i in 2:md[d]){
      Omegad[[d]][i,i-1]=-theta
      Omegad[[d]][i-1,i]=-theta
    }
    diag(Omegad[[d]]) <- 1+theta^2
    
    
    
    #Derivada
    
    OmegadPrima[[d]] <- matrix(0,nrow=md[d],ncol=md[d])
    for(i in 2:md[d]){
      OmegadPrima[[d]][i,i-1]=-1
      OmegadPrima[[d]][i-1,i]=-1
    }
    diag(OmegadPrima[[d]]) <- 2*theta
    
    
    
    Xd[[d]] <- X[a[[d]],]
    Ved <- diag(sigma2edi[a[[d]]])                                           # Matriz Ved
    Vd <- (sigmau1 *(uno%*%t(uno)) + sigmau2 * Omegad[[d]] +  Ved)                # Matriz de varianza
    Vd.inv[[d]] <- solve(Vd)                                                 # Inversa de la matriz de varianza en d submatrices
    Ved.inv[[d]] <- solve(Ved)                                               # Inversa de la matriz Sigma_ed en d submatrices 
    VinvOmega[[d]] <-  Vd.inv[[d]]%*%  Omegad[[d]]                           # Producto de V^-1_d por Omega 
    OmegaVinvOmega[[d]] <- t(VinvOmega[[d]])%*%Omegad[[d]]                   # Producto de Omega por V^-1_d por Omega
    VedinvXd[[d]] <- Ved.inv[[d]]%*%Xd[[d]]                                  # Producto de V^-1_ed  por  X_d para las d submatrices
    VinvOmegaPrima[[d]] <-  Vd.inv[[d]]%*%OmegadPrima[[d]]                   # Producto de V^-1_d por OmegaPrima 
    OmegaVinvOmegaPrima[[d]] <- Omegad[[d]]%*%Vd.inv[[d]]%*%OmegadPrima[[d]] # Producto de Omega por V^-1_d por OmegaPrima 
    OmegaPrimaVinvOmegaPrima[[d]] <- OmegadPrima[[d]]%*%Vd.inv[[d]]%*%OmegadPrima[[d]]   # Producto de OmegaPrima por V^-1_d por OmegaPrima 
    
    
    T11d[[d]] <- sigmau1-sigmau1^2*t(uno)%*%Vd.inv[[d]]%*%uno
    T12d[[d]] <- -sigmau1*sigmau2*t(uno)%*%VinvOmega[[d]]
    T22d[[d]] <- sigmau2*Omegad[[d]] -sigmau2^2*OmegaVinvOmega[[d]] 
    
    #OmegaVinvOmegadSinvXd[[d]] <- OmegaVinvOmega[[d]]%*%SinvXd[[d]]          # Producto Omegad por V^-1_d  por  Omegad por Sigma^-1_ed  por  X_d para las d submatrices       
    
    
    g1.a[[d]] <- uno%*%T11d[[d]]%*%t(uno) + uno%*%T12d[[d]] + t(T12d[[d]])%*%t(uno) + T22d[[d]]
    
    g2.a[[d]] <- Xd[[d]] - (uno%*%T11d[[d]]%*%t(uno) + uno%*%T12d[[d]] + t(T12d[[d]])%*%t(uno) + T22d[[d]])%*% VedinvXd[[d]]   # Primera parte de g2 (la segunda es su traspuesta)
    
    
    auxq[[d]] <- sigmau1*unotuno + sigmau2* Omegad[[d]]
    unotunoVinv[[d]] <- unotuno%*%Vd.inv[[d]]
    
    q11[[d]] <- unotunoVinv[[d]]%*%unotuno - 2*unotunoVinv[[d]]%*%unotunoVinv[[d]]%*%auxq[[d]]+auxq[[d]]%*%Vd.inv[[d]]%*%unotunoVinv[[d]]%*%unotunoVinv[[d]]%*%auxq[[d]]
    
    
    q12[[d]]  <- unotuno%*%VinvOmega[[d]]-unotunoVinv[[d]]%*%t(VinvOmega[[d]])%*%auxq[[d]] - auxq[[d]]%*%t(unotunoVinv[[d]])%*%VinvOmega[[d]] + auxq[[d]]%*%t(unotunoVinv[[d]])%*%VinvOmega[[d]]%*%Vd.inv[[d]]%*%auxq[[d]]
    
    q13[[d]] <- sigmau2*unotuno%*%VinvOmegaPrima[[d]] - sigmau2*unotunoVinv[[d]]%*%t(VinvOmegaPrima[[d]])%*%auxq[[d]] - sigmau2*auxq[[d]]%*%t(unotunoVinv[[d]])%*%VinvOmegaPrima[[d]] 
    + sigmau2*auxq[[d]]%*%Vd.inv[[d]]%*%unotunoVinv[[d]]%*%t(VinvOmegaPrima[[d]]%*%auxq[[d]])
    
    
    q22[[d]]  <- OmegaVinvOmega[[d]] - 2*OmegaVinvOmega[[d]]%*%Vd.inv[[d]]%*%auxq[[d]] + auxq[[d]]%*%Vd.inv[[d]]%*%OmegaVinvOmega[[d]]%*%Vd.inv[[d]]%*%auxq[[d]]    
    
    q23[[d]] <- sigmau2*OmegaVinvOmegaPrima[[d]]  - sigmau2*OmegaVinvOmegaPrima[[d]] %*%Vd.inv[[d]]%*%auxq[[d]] - sigmau2*auxq[[d]]%*%VinvOmega[[d]]%*%VinvOmegaPrima[[d]] 
    + sigmau2*auxq[[d]]%*%Vd.inv[[d]]%*%OmegaVinvOmegaPrima[[d]]%*%Vd.inv[[d]]%*%auxq[[d]]	
    q33[[d]] <- sigmau2^2*OmegaPrimaVinvOmegaPrima[[d]] - 2*sigmau2^2*OmegaPrimaVinvOmegaPrima[[d]]%*%Vd.inv[[d]]%*%auxq[[d]] + sigmau2^2*auxq[[d]]%*%Vd.inv[[d]]%*%OmegaPrimaVinvOmegaPrima[[d]]%*%Vd.inv[[d]]%*%auxq[[d]]
    
    
    Q.inv <- Q.inv + t(Xd[[d]])%*%Vd.inv[[d]]%*%Xd[[d]]                      # Inversa de Q. Posteriormente se calcula Q
  }
  
  Q <- solve(Q.inv)
  
  
  # Calgulo de g
  
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
