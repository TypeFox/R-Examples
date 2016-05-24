mse.saery.indep <-
function(X, D, md, sigma2edi, sigmau1, sigmau2, Fsig) {
  p <- ncol(X)
  a <- list(1:md[1])
  mdcum <- cumsum(md)
  for(d in 2:D)
    a[[d]] <- (mdcum[d-1]+1):mdcum[d]
  Xd <- Vd.inv <- Ved.inv <- T11d <- T12d <- T22d <- VedinvXd <- auxq <- V1dVinv <- g1.a <- g2.a <- q11 <- q12 <- q22 <- list()
  Q.inv <- matrix(0, nrow=p, ncol=p)
  
  for(d in 1:D) {
    uno <- rep(1,md[d])
    V1d <- uno%*%t(uno)
    V2d <- diag(uno)
    Xd[[d]] <- X[a[[d]],]
    Ved <- diag(sigma2edi[a[[d]]])
    Vd <- (sigmau1*V1d + sigmau2*V2d + Ved)
    Vd.inv[[d]] <- solve(Vd)
    Ved.inv[[d]] <- solve(Ved)
    VedinvXd[[d]] <- Ved.inv[[d]]%*%Xd[[d]]
    
    T11d[[d]] <- sigmau1-sigmau1^2*t(uno)%*%Vd.inv[[d]]%*%uno
    T12d[[d]] <- -sigmau1*sigmau2*t(uno)%*%Vd.inv[[d]]
    T22d[[d]] <- diag(sigmau2,md[d],md[d])-sigmau2^2*Vd.inv[[d]]
    
    g1.a[[d]] <- uno%*%T11d[[d]]%*%t(uno) + uno%*%T12d[[d]] + t(T12d[[d]])%*%t(uno) + T22d[[d]]
    g2.a[[d]] <- Xd[[d]] - (uno%*%T11d[[d]]%*%t(uno) + uno%*%T12d[[d]] + t(T12d[[d]])%*%t(uno) + T22d[[d]])%*% VedinvXd[[d]]
    auxq[[d]] <- sigmau1*V1d + sigmau2*V2d
    V1dVinv[[d]] <- V1d%*%Vd.inv[[d]]    
    q11[[d]] <- V1dVinv[[d]]%*%V1d - 2*V1dVinv[[d]]%*%V1dVinv[[d]]%*%auxq[[d]] + auxq[[d]]%*%Vd.inv[[d]]%*%V1dVinv[[d]]%*%V1dVinv[[d]]%*%auxq[[d]]
    q12[[d]] <- V1dVinv[[d]] - V1dVinv[[d]]%*%Vd.inv[[d]]%*%auxq[[d]] - auxq[[d]]%*%Vd.inv[[d]]%*%V1dVinv[[d]] + auxq[[d]]%*%Vd.inv[[d]]%*%V1dVinv[[d]]%*%Vd.inv[[d]]%*%auxq[[d]]
    q22[[d]] <- Vd.inv[[d]] - 2*Vd.inv[[d]]%*%Vd.inv[[d]]%*%auxq[[d]] + auxq[[d]]%*%Vd.inv[[d]]%*%Vd.inv[[d]]%*%Vd.inv[[d]]%*%auxq[[d]]
    Q.inv <- Q.inv + t(Xd[[d]])%*%Vd.inv[[d]]%*%Xd[[d]]
  }
  Q <- solve(Q.inv)
  
  g1 <- g2 <- g3 <- list()
  for(d in 1:D){
    g1[[d]] <- diag(g1.a[[d]])
    g2[[d]] <- diag(g2.a[[d]]%*%Q%*%t(g2.a[[d]]))
    q11[[d]] <- diag(q11[[d]])
    q12[[d]] <- diag(q12[[d]])
    q22[[d]] <- diag(q22[[d]])
  }
  for(d in 1:D){
    g3[[d]] <- vector()
    for(i in 1:md[d]){
      g3[[d]][i] <- sum(diag(matrix(c(q11[[d]][i],q12[[d]][i],q12[[d]][i],q22[[d]][i]),nrow=2)%*%solve(Fsig)))
    }
  }   
  
  g1 <- unlist(g1)
  g2 <- unlist(g2)
  g3 <- unlist(g3)
  
  return(g1 + g2 + 2*g3)
  
}
