simrel <-
function(n, p, m, q, relpos, gamma, R2, ntest=NULL, muY=NULL, muX=NULL,sim=NULL){
  
  if(q<m) stop(paste("the number of relevant predictors must at least be equal to",m,"\n"))
  if(!(is.null(sim))){
    if((sim$m!=m)||(sim$p!=p)||(sim$relpos!=relpos)||(sim$q!=q)||(sim$gamma!=gamma)||(sim$R2!=R2)){
      stop("Different parameter setting than previous object")
    }
  }
  if(length(relpos)!=m){stop("Mismatch between number of relevant components and the given positions in 'relpos'\n")}
  
  irrelpos <- (1:p)[-relpos]
  extradim <- q - m
  if(is.null(sim)){
    qpos <- sort(c(relpos, sample(irrelpos, extradim, replace=F)))
  }else{
    qpos <- sim$relpred
    }
  lambdas <- exp(-gamma*(1:p))/exp(-gamma)
  
  #Construction of Sigma
  SigmaZ <- diag(lambdas)   
  
  SigmaZinv <- diag(1/lambdas)
  Sigmazy <- matrix(0,p,1)
  if(is.null(sim)){
    r <- runif(m, 0, 1)*sample(c(-1,1),m,replace=TRUE)
  }else{
    r <- sim$r
  }
  Sigmazy[relpos,] <- sign(r)*sqrt(R2*abs(r)/sum(abs(r))*lambdas[relpos])
  SigmaY <- 1
  Sigma <- rbind(c(SigmaY,t(Sigmazy)),cbind(Sigmazy, SigmaZ))
  
  #Checking that Sigma is PD
  pd <- all(eigen(Sigma)$values>0)
  
  #Finding a rptation matrix which rotates the latent components to yield
  #relevant X-predictors
  if(is.null(sim)){
    Q <- matrix(rnorm(q^2),q)
    Q <- scale(Q,scale=F)
    qrobj <- .QR(Q)
    Rq <- qrobj$Q
    R <- diag(p)
    R[qpos, qpos] <- Rq
    if(q<(p-1)){  
      Q <- matrix(rnorm((p-q)^2),(p-q))
      Q <- scale(Q,scale=F)
      qrobj <- .QR(Q)
      Rnq <- qrobj$Q
      R[(1:p)[-qpos],(1:p)[-qpos]] <- Rnq
    }
  }else{
    R <- sim$Rotation
  }
  
  #The true regression coefficients
  betaZ <- SigmaZinv%*%Sigmazy
  betaX <- R%*%betaZ
  beta0 <- 0
  if(!(is.null(muY))){beta0 <- beta0 + muY}
  if(!(is.null(muX))){beta0 <- beta0 - t(betaX)%*%muX}
  #The (true) coefficient of determination, R^2
  R2 <- t(Sigmazy)%*%betaZ
  #Minimum prediction error
  minerror <- SigmaY - R2
  
  
  #Simulating training and test data
  if(pd){
    Sigmarot<-chol(Sigma)
    Ucal<-matrix(rnorm(n*(p+1),0,1),nrow=n)
    U1cal<-Ucal%*%Sigmarot
    Y <- U1cal[,1,drop=F]
    if(!(is.null(muY))){Y <- Y + rep(muY,n)}
    Z <- U1cal[,2:(p+1)]
    X <- Z%*%t(R)
    if(!(is.null(muX))){X <- sweep(X, 2, muX, "+")}
    colnames(X) <- as.character(1:p)
    
    #Testdata
    if(!is.null(ntest)){
      Utest<-matrix(rnorm(ntest*(p+1),0,1),nrow=ntest)
      U1test<-Utest%*%Sigmarot
      TESTY<-U1test[,1,drop=F]
      if(!(is.null(muY))) TESTY <- TESTY + rep(muY,ntest)
      TESTZ <- U1test[,2:(p+1)]
      TESTX <- TESTZ%*%t(R)
      if(!(is.null(muX))){TESTX <- sweep(TESTX, 2, muX, "+")}
      colnames(TESTX) <- as.character(1:p)
    }else{
      TESTX <- NULL
      TESTY <- NULL
    }
    
   
  }else{stop("Correlation matrix is not positive definit \n")}
  
  res <- list()
  res$call <- match.call()
  res$X <- X
  res$Y <- Y
  res$beta <- betaX
  res$beta0 <- beta0
  res$relpred <- qpos
  res$TESTX <- TESTX
  res$TESTY <- TESTY
  res$n <- n
  res$p <- p
  res$m <- m
  res$q <- q
  res$gamma <- gamma
  res$lambda <- lambdas
  res$R2 <- drop(R2)
  res$relpos <- relpos
  res$minerror <- minerror
  res$r <- r
  res$Sigma <- Sigma
  res$Rotation <- R
  class(res) <- "simrel"
  return(res)
}
