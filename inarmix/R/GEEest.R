GEEests <- function(weights,pars,XX,yy,nvec,R.alpha,A1,A2,E1,E2,E3,E4)  {
  beta <- as.numeric(pars$beta)
  alpha <- pars$alpha
  gamma <- pars$gamma
  eps <- 1e-4
  count <- 0
  #### Iterate GEEs at most 10 times
  GEEiter <- 10
  for(k in 1:GEEiter)   {
    count <- count + 1
    R.alpha <- (alpha*E1 + (1 + alpha^2)*E2 + E3)/(1 - alpha^2) + E4 
    
    newbeta <- UpdateBeta(betinit=beta, XX=XX,yy=yy,A1.store=A1,A2.store=A2,R.alpha=R.alpha,nvec,weights)
    gamma <- UpdateGamma(XX,yy,A1,R.alpha,newbeta,nvec,weights)
    newalpha <- UpdateAlpha(XX,yy,A.store=A1,R.alpha=R.alpha,newbeta,alpha,gamma,nvec=nvec,E1=E1,E2=E2,E3=E3,weights) 
    
    done <- (max(abs(newbeta - beta)/(abs(newbeta) + 1)) < eps)
    if (done)  {
      beta <- newbeta
      alpha <- max(alpha,.01)
      gamma <- max(gamma,.01)
      break
    }
    
    beta <- newbeta   
    alpha <- newalpha 
  }
  conv <- ifelse(count<GEEiter,1,0)
  return(list(beta=newbeta,alpha=newalpha,gamma=gamma,conv = conv))
}


UpdateGamma <- function(XX,yy,A.store,R.alpha,beta,nvec,weights)  {
  ### N.weight is the N_{i} in the notes
  expanded.weights <- rep(weights,nvec)
  mu.stack <- as.vector(exp(XX%*%beta))
  diag(A.store) <- (1/sqrt(mu.stack))     
  weighted.dat <- sqrt(expanded.weights)*(yy - mu.stack)
  N.weight <- sum(weights*nvec)
  
  gam <- (1/N.weight)*crossprod(A.store%*%weighted.dat,R.alpha%*%A.store%*%weighted.dat) - 1
  
  if(as.numeric(gam) <= 0)  {
    print("The scale parameter is near the boundary")
    gam <- .01
  }
  
  return(as.numeric(gam))
}


UpdateAlpha <- function(XX,yy,A.store,R.alpha,beta,alpha,gamma,nvec,E1,E2,E3,weights)  {
  expanded.weights <- rep(weights,nvec)
  mu.stack <- as.vector(exp(XX%*%beta))
  diag(A.store) <- (1/sqrt(mu.stack))
  weighted.dat <- sqrt(expanded.weights)*(yy - mu.stack)
  N.weight <- sum(weights*(nvec-1))
  
  eps <- .001
  ### Perform a maximum of 10 Fisher scoring iterations 
  for(k in 1:20)  {
    R.alpha <-  -((1 + alpha^2)*E1 + 4*alpha*E2 + 2*alpha*E3)/((1 - alpha^2)^2) 
    
    LHT <- as.numeric(crossprod(A.store%*%weighted.dat,R.alpha%*%A.store%*%weighted.dat))
    RHT <- (2*alpha*N.weight*(1+gamma))/(1 - alpha^2)
    
    alpha.der <- (2*(1 + alpha^2)*N.weight*(1+gamma))/((1 - alpha^2)^2)
    
    newalpha <- alpha + (LHT + RHT)/alpha.der
    
    if((newalpha < 0) | (newalpha >= 1)) {
      print("The autocorrelation parameter is near the boundary")
      if(newalpha < 0) {
        newalpha <- .01
      }
      else {
        newalpha <- .95
      }
    }
    
    done <- (abs(newalpha - alpha)/(abs(newalpha) + 1) < eps)
    if (done)  {
      alpha <- newalpha
      break
    }
    alpha <- newalpha
  }
  return(newalpha)
}

UpdateBeta <- function(betinit,XX,yy,A1.store,A2.store,R.alpha,nvec,weights)  {
  expanded.weights <- rep(weights,nvec)
  beta <- betinit
  eps <- .01
  #### perform a maximum of 30 Fisher scoring iterations
  for (k in 1:30)  {
    mu.stack <- as.vector(exp(XX%*%beta))
    diag(A1.store) <- sqrt(expanded.weights*mu.stack)
    diag(A2.store) <- (1/sqrt(mu.stack))
    weighted.dat <- sqrt(expanded.weights)*(yy - mu.stack)
    
    beta.esteq <- as.vector(crossprod(A1.store%*%XX,R.alpha%*%A2.store%*%weighted.dat))
    beta.hess <- crossprod(A1.store%*%XX,R.alpha%*%A1.store%*%XX)
    
    if(rcond(beta.hess) < sqrt(.Machine$double.eps)) {
      beta.hess <- beta.hess + 2e-8*diag(nrow(beta.hess)) 
    }
    newbeta <- beta + solve(beta.hess,beta.esteq)
    
    done <- (max(abs(newbeta - beta)/(abs(newbeta) + 1)) < eps)
    if (done)  {
      beta <- newbeta
      break
    }
    beta <- newbeta
  }
  return(newbeta)
}


BuildBasis <- function(nvec)  {
  N <- sum(nvec)
  m <- length(nvec)
  a1 <- a2 <- a3 <- a4 <- vector("numeric", N)
  index <- c(cumsum(nvec) - nvec, N)
  for (i in 1:m)  {
    if(nvec[i] > 1)  {
      a1[(index[i]+1) : index[i+1]] <- c(rep(-1,times = nvec[i]-1),0)
      a2[(index[i]+1) : index[i+1]] <- c(0,rep(1,times=nvec[i]-2),0)
      a3[(index[i]+1) : index[i+1]] <- c(1,rep(0,times=nvec[i]-2),1)
      a4[(index[i]+1) : index[i+1]] <- c(rep(0,times=nvec[i]))
    }
    else if (nvec[i] == 1)  {
      a1[(index[i]+1) : index[i+1]] <- 0
      a2[(index[i]+1) : index[i+1]] <- 0
      a3[(index[i]+1) : index[i+1]] <- 0
      a4[(index[i]+1) : index[i+1]] <- 1
    }
  }
  
  ar.pieces <- list(a1, rep(0,N))
  E1 <- bandSparse(N,k=c(-1,0),diagonals = c(ar.pieces),symmetric=TRUE)
  E2 <- Diagonal(x = a2)
  E3 <- Diagonal(x = a3)
  E4 <- Diagonal(x = a4)
  
  return(list(E1=E1,E2=E2,E3=E3,E4=E4))
}



ConstructIndSet <- function(z,nclasses,ncum)  {
  #### z is a vector indicating class membership (e.g. z = (1,2,1,4,3,3))
  
  #### This function returns a matrix whose (i,j) element is given by
  #### A[i,j] = 1 if the observation j belongs to class i
  
  #### For example, if there are two subjects each with two observations
  #### and the first subject belongs to class 1 and the second subject belongs to class 2, then
  #### indset = rbind(c(1,1,0,0),c(0,0,1,1))
  
  m <- length(ncum) - 1
  N <- max(ncum) - 1
  indset <- matrix(0,nrow=nclasses,ncol=N)
  
  for(k in 1:nclasses)  {
    for (j in 1:m)  {
      ind1 <- ncum[j]
      ind2 <- ncum[j+1] - 1
      val <- ifelse(z[j]==k,1,0)
      indset[k,ind1:ind2] <- val
    }
  }
  return(indset)
}