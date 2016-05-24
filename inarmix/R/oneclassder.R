OuterProd <- function(beta,alpha,gamma,nvec,ylist,Xlist)  {
  m <- length(nvec)
  p <- length(beta) + 2
  
  const1 <- (-2*alpha)/(1 - alpha^2)
  const2 <- (2*(1 + alpha^2))/((1 - alpha^2)^2)
  
  outer.list <- list()
  DD <- matrix(0,nrow=p,ncol=p)
  stack.mat <- der.mat <- matrix(0,nrow=p^2,ncol=m)
  for(i in 1:m)  {
    X.tmp <- Xlist[[i]]
    
    yy <- ylist[[i]]
    Rtmp <- R.compute(alpha,T=length(yy))
    R.inv <- Rtmp$R.inv
    R.sand <- Rtmp$R.sand
    T <- length(yy)
    
    mu.vec <- exp(X.tmp%*%beta)
    tt <- as.vector(1/sqrt(mu.vec))
    ttt <- as.vector(sqrt(mu.vec))
    
    A <- diag(tt)
    A2 <- diag(ttt)
    
    u1 <- crossprod(A2%*%X.tmp,R.inv%*%A%*%(yy - mu.vec))/(1 + gamma)
    u2 <- crossprod(A%*%(yy - mu.vec),R.sand%*%A%*%(yy - mu.vec)) + (2*(1 + gamma)*alpha*(T-1))/(1 - alpha^2)
    u2 <- u2/(1 + gamma)
    u3 <- crossprod(A%*%(yy - mu.vec),R.inv%*%A%*%(yy - mu.vec)) - T*(1 + gamma)
    u3 <- u3/((1 + gamma)^2) 
    
    a <- c(u1,u2,u3)
    outer.list[[i]] <- outer(a,a)
    
    stack.mat[,i] <- c(outer(a,a))
    
    times <- c(1:T)
    H <- abs(outer(times, times, "-"))
    R <- alpha^H
    
    len.beta <- length(beta)
    A.1 <- crossprod(A2%*%X.tmp,R.inv%*%A2%*%X.tmp)/(1 + gamma) 
    A.2 <- matrix(0,nrow=2,ncol=len.beta)
    A.3 <- matrix(0,nrow=2,ncol=2)
    
    for(k in 1:len.beta)  {
      X.special <- outer(X.tmp[,k],X.tmp[,k],FUN="+")/2                
      Var.mat <- (A2%*%R%*%A2)
      V1 <- X.special*(Var.mat)
      inv.alpha <- (A%*%R.sand%*%A)
      inv.gamma <- (A%*%R.inv%*%A)/(1 + gamma)
      
      A.2[1,k] <- sum(V1*inv.alpha) 
      A.2[2,k] <- sum(V1*inv.gamma)
    }
    A.3[1,1] <- const2*(T-1)
    A.3[1,2] <- A.3[2,1] <- (const1*(T-1))/(1 + gamma)
    A.3[2,2] <- T/((1 + gamma)^2) 
    A.4 <- matrix(0,nrow=len.beta,ncol=2)
    
    Der.mat <- cbind(rbind(A.1,A.2),rbind(A.4,A.3))
    DD <- DD + Der.mat  
  }
  DD <- DD/m
  
  tmp <- apply(stack.mat,1,mean)
  FinMat <- matrix(tmp,nrow=p,ncol=p)/m
  
  res <- list()
  res$V <- FinMat
  res$D <- DD
  return(res)
}


R.compute <- function(alpha,T)  {
  ### This function computes the inverse correlation matrix
  ### for a given autocorrelation value and observational length.
  ### It also computes the "R sandwich" matrix.
  ### R.sand = R.inv*dR*R.inv = -dR.inv/dalpha.
  
  if(T > 2)  {
    det1 <- 1/(1 - alpha^2)
    det2 <- -(det1^2)    
    
    R1 <- R2 <- matrix(0,nrow=T,ncol=T)
    
    R1[1,1] <- R1[T,T] <- 1
    R1[1,2] <- R1[T,T-1] <- -alpha
    
    R2[1,1] <- R2[T,T] <- 2*alpha
    R2[1,2] <- R2[T,T-1] <- -(1 + alpha^2)
    
    for(i in 2:(T-1))  {
      R1[i,i-1] <- R1[i,i+1] <- -alpha
      R1[i,i] <- 1 + alpha^2
      
      R2[i,i-1] <- R2[i,i+1] <- -(1 + alpha^2)
      R2[i,i] <- 4*alpha
    }
    
    res <- list()
    res$R.inv <- det1*R1
    res$R.sand <- det2*R2
    return(res) 
  }
  else if (T==2)  {
    R1 <- R2 <- matrix(0,2,2)
    R1[1,1] <- R1[2,2] <- 1/(1 - alpha^2)
    R1[1,2] <- R1[2,1] <- -alpha/(1 - alpha^2)
    
    R2[1,1] <- R2[2,2] <- (-2*alpha)/((1 - alpha^2)^2)
    R2[1,2] <- R2[2,1] <- (1 + alpha^2)/((1 - alpha^2)^2)
    
    res <- list() 
    res$R.inv <- R1
    res$R.sand <- R2
    return(res)
  }
  else if(T==1)  {
    res <- list()
    res$R.inv <- 1
    res$R.sand <- 0
    return(res)
  }
}
