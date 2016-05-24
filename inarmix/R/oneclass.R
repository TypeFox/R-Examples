OneclassEsts <- function(yy,XX,nvec,m,initvals=NULL)  {
  #### Initialize parameter values
  if(is.null(initvals))  {
    init.intercept <- log(max(mean(yy),.01))
    if(ncol(XX)==1)   {
      beta <- init.intercept
    }
    else  {
      beta <- c(init.intercept,rep(0,times=ncol(XX) - 1))
    }
    alpha <- .2
  }
  else  {
    beta <- initvals$beta
    alpha <- initvals$alpha
    gamma <- initvals$gamma
  }
  
  #### set up Matrix storage
  N <- sum(nvec)
  A1.st <- Diagonal(N,x=rep(1,N))
  A2.st <- Diagonal(N,x=rep(1,N))
  A.st <- Diagonal(N,x=rep(1,N))
  #### builds it almost immediately
  
  #### set up R.alpha "basis" matrices
  tmp <- BuildBasis_oneclass(nvec)
  E1.st <- tmp$E1
  E2.st <- tmp$E2
  E3.st <- tmp$E3
  E4.st <- tmp$E4
  
  ### iterate estimating equations a maximum of 20 times.
  eps <- 1e-4
  count <- 1
  for(k in 1:20)   {
    R.alpha <- (alpha*E1.st + (1 + alpha^2)*E2.st + E3.st)/(1 - alpha^2) + E4.st 
    ### pretty quick once stored (about .02 seconds)
    
    
    newbeta <- UpdateBeta_oneclass(betinit=beta, XX=XX,yy=yy,A1.store=A1.st,A2.store=A2.st,R.alpha=R.alpha)
    gamma <- UpdateGamma_oneclass(XX,yy,A.st,R.alpha,newbeta)
    newalpha <- UpdateAlpha_oneclass(XX,yy,A.store=A.st,R.alpha=R.alpha,newbeta,alpha,gamma,nvec=nvec,E1=E1.st,E2=E2.st,E3=E3.st) 
    
    
    done <- (max(abs(newbeta - beta)/(abs(newbeta) + 1)) < eps)
    if (done)  {
      beta <- newbeta   
      alpha <- newalpha 
      break
    }
    beta <- newbeta   
    alpha <- newalpha 
  }
  
  results <- list() 
  results$beta <- as.vector(newbeta)
  results$alpha <- alpha
  results$gamma <- gamma
  return(results)    
}


UpdateBeta_oneclass <- function(betinit,XX,yy,A1.store,A2.store,R.alpha)  {
  beta <- betinit
  eps <- .01
  #### perform a maximum of 20 Fisher scoring iterations
  for (k in 1:20)  {
    mu.stack <- as.vector(exp(XX%*%beta))
    diag(A1.store) <- sqrt(mu.stack)
    diag(A2.store) <- 1/sqrt(mu.stack)
    
    beta.esteq <- as.vector(crossprod(A1.store%*%XX,R.alpha%*%A2.store%*%(yy - mu.stack)))
    beta.hess <- crossprod(A1.store%*%XX,R.alpha%*%A1.store%*%XX)
    
    if(rcond(beta.hess) < sqrt(.Machine$double.eps)) {
      beta.hess <- beta.hess + 2e-8*diag(nrow(beta.hess)) 
    }
    newbeta <- beta + solve(beta.hess,beta.esteq)
    
    done <- (max(abs(newbeta - beta)/(abs(newbeta) + 1)) < eps)
    if (done)  {break}
    beta <- newbeta
  }
  return(newbeta)
}

UpdateGamma_oneclass <- function(XX,yy,A.store,R.alpha,beta)  {
  mu.stack <- as.vector(exp(XX%*%beta))
  diag(A.store) <- 1/sqrt(mu.stack)     
  nn <- length(yy)
  
  gam <- (1/nn)*crossprod(A.store%*%(yy - mu.stack),R.alpha%*%A.store%*%(yy - mu.stack)) - 1
  return(as.numeric(gam))
}

UpdateAlpha_oneclass <- function(XX,yy,A.store,R.alpha,beta,alpha,gamma,nvec,E1,E2,E3)  {
  m <- length(nvec)
  N <- sum(nvec)
  mu.stack <- as.vector(exp(XX%*%beta))
  diag(A.store) <- 1/sqrt(mu.stack)
  
  eps <- .001
  ### Perform a maximum of 10 Fisher scoring iterations 
  for(k in 1:10)  {
    R.alpha <-  -((1 + alpha^2)*E1 + 4*alpha*E2 + 2*alpha*E3)/((1 - alpha^2)^2) 
    
    LHT <- as.numeric(crossprod(A.store%*%(yy - mu.stack),R.alpha%*%A.store%*%(yy - mu.stack)))
    RHT <- (2*alpha*(N-m)*(1+gamma))/(1 - alpha^2)
    
    alpha.der <- (2*(1 + alpha^2)*(N-m)*(1+gamma))/((1 - alpha^2)^2)
    
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
    if (done)  {break}
    alpha <- newalpha
  }
  return(newalpha)
}


BuildBasis_oneclass <- function(nvec)  {
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