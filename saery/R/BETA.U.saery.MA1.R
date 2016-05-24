BETA.U.saery.MA1 <-
function(X, ydi, D, md, sigma2edi, sigmau1, sigmau2, theta) {
  
  p <- ncol(X)
  a <- list(1:md[1])
  mdcum <- cumsum(md)
  for(d in 2:D)
    a[[d]] <- (mdcum[d-1]+1):mdcum[d]
  
  yd <- Xd <- Vd <- Vd.inv <- list()
  Q.inv <- matrix(0, nrow=p, ncol=p)
  XVy <- 0
  for(d in 1:D) {
    uno<-rep(1,md[d])
    yd[[d]] <- ydi[a[[d]]]
    Xd[[d]] <- X[a[[d]],]
    
    Omegad <- matrix(0,nrow=md[d],ncol=md[d])
    for(i in 2:md[d]){
      Omegad[i,i-1]<- -theta
      Omegad[i-1,i]<- -theta
    }
    diag(Omegad) <- 1+theta^2
    
    
    
    
    Vd <- (sigmau1 *(uno%*%t(uno)) + sigmau2 * Omegad + diag(sigma2edi[a[[d]]]))          
    Vd.inv[[d]] <- solve(Vd)                              
    Q.inv <- Q.inv + t(Xd[[d]])%*%Vd.inv[[d]]%*%Xd[[d]]         
    XVy <- XVy + t(Xd[[d]])%*%Vd.inv[[d]]%*%yd[[d]]             
  }
  Q <- solve(Q.inv)
  
  beta <- Q%*%XVy
  
  
  u1 <- u2 <- list()
  for(d in 1:D){ 
    u1[[d]] <- sigmau1*(t(uno)%*%Vd.inv[[d]]%*%(yd[[d]]-Xd[[d]]%*%beta))
    u2[[d]] <- sigmau2*(Omegad%*%Vd.inv[[d]]%*%(yd[[d]]-Xd[[d]]%*%beta))
  }
  
  u1 <- as.matrix(unlist(u1))
  u2 <- as.matrix(unlist(u2))
  
  
  
  return(list(beta,u1,u2))
}
