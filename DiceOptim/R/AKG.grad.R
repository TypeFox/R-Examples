AKG.grad <- function(x, model, new.noise.var=0, type = "UK", envir=NULL){

  ########## Convert x in proper format(s) ###
  d <- length(x)
  if (d != model@d){ stop("x does not have the right size") }
  newdata.num <- as.numeric(x)
  newdata <- data.frame(t(newdata.num))
  colnames(newdata) = colnames(model@X)
  tau2.new <- new.noise.var
 
  ######### Intermediate values ##########
  T <- model@T
  z <- model@z
  U <- model@M
  F.X <- model@F
  
  if (is.null(envir))
  {
    ######### Compute prediction at x #########
    predx <- predict(model, newdata=newdata, type=type, checkNames = FALSE)
    mk.x <- predx$mean
    sk.x <- predx$sd
    c.x  <- predx$c
    V.x <- predx$Tinv.c
    F.x <- model.matrix(model@trend.formula, data=newdata)
    
    if (sk.x > sqrt(model@covariance@sd2)/1e6 || model@covariance@sd2 > 1e-20)
    {
      ######### Compute prediction at X #########
      predX <- predict.km(model, newdata=model@X, type=type, checkNames = FALSE)
      mk.X <- predX$mean
      V.X <- predX$Tinv.c
      
      ######### Compute cn #########
      if (type=="UK")
      { tuuinv <- solve(t(U)%*%U)
        mu.x <- (F.X - t(V.X)%*%U)%*% tuuinv %*% t(F.x - t(V.x)%*%U)
      } else {tuuinb <- mu.x <- 0}
      cn <- c.x - t(V.X)%*%V.x + mu.x
      cn <- c(cn, sk.x^2)
      
      ######### Compute a and b #########
      a <- c(mk.X, mk.x)
      b <- cn / sqrt(tau2.new + sk.x^2)
      sQ <- b[model@n+1]
    }
  } else
  { 
    ######## Get all necessary data #######
    toget <- matrix(c("mk.x", "c.x", "V.x", "sk.x", "F.x", "tuuinv", "mk.X", "V.X", "mu.x", "cn", "sQ", 
                      "Isort", "Iremove", "A1", "at", "bt", "ct"),1,17)
    apply(toget, 2, get, envir=envir)
    mk.x   <- envir$mk.x  #predx$sd
    c.x    <- envir$c.x   #predx$c
    V.x    <- envir$V.x   #predx$Tinv.c
    sk.x   <- envir$sk.x  #predx$sd
    F.x    <- envir$F.x
    tuuinv <- envir$tuuinv
    mk.X   <- envir$mk.X
    V.X    <- envir$V.X   #predX$Tinv.c
    mu.x   <- envir$mu.x  
    cn     <- envir$cn
    sQ     <- envir$sQ    #b[model@n+1]
  }
  if (sk.x < sqrt(model@covariance@sd2)/1e6 || model@covariance@sd2 < 1e-20)
  { grad.AKG <- rep(0,model@d)
  }  else
  {
    ######### Some derivatives ########## 
    dc.x <- covVector.dx(model@covariance, x=newdata.num, X=model@X, c=c.x)
    f.deltax <- trend.deltax(x=newdata.num, model=model)
    W <- backsolve(t(T), dc.x, upper.tri=FALSE)
    
    ######### Gradient of a ##########
    a.grad <- matrix(0,model@d,model@n+1)
    a.grad[,model@n+1] <- t(W)%*%z + t(model@trend.coef%*%f.deltax)

    ######### Gradient of b ##########
    b.grad <- matrix(0,model@d, model@n+1)
    
    if (type=="UK")
    { sk2.x.grad <-  t( -2*t(V.x)%*%W + 2*(F.x - t(V.x)%*%U )%*% tuuinv %*%
                        (f.deltax - t(t(W)%*%U) ))    
      mu.x.grad  <- t(2*(F.X - t(V.X)%*%U )%*% tuuinv %*% (f.deltax - t(t(W)%*%U) ))
    } else
    { sk2.x.grad <- t( -2*t(V.x)%*%W)
      mu.x.grad  <- matrix(0,model@d, model@n)
    }
#     print(t(dc.x) - t(t(V.X)%*%W))
    cn.grad <- t(dc.x) - t(t(V.X)%*%W) + mu.x.grad
    
    b.grad[,1:model@n] <- cn.grad/sqrt(tau2.new+sk.x^2) - sk2.x.grad%*%cn[1:model@n]*as.numeric(1/2/(tau2.new+sk.x^2)^(3/2))
    b.grad[,model@n+1] <- as.numeric(sk.x^2*(2*tau2.new+sk.x^2)/(tau2.new+sk.x^2)^2)*sk2.x.grad/as.numeric(2*sQ)
    
    ######### Careful: the AKG is written for MAXIMIZATION #########
    ######### Minus signs have been added where necessary ##########
    
    if (is.null(envir))
    {
      #--------------- Recompute at, bt, ct ----------------------------
      A <- -a
      B <- b
      B.grad <- b.grad
      A.grad <- -a.grad
      
      ######## Sort and reduce A and B ####################
      # Sort by increasing B values
      Isort <- order(x=B,y=A)
      b <- B[Isort]	
      b.grad <- B.grad[,Isort]
      a <- A[Isort]	
      a.grad <- A.grad[,Isort]
      
      # Remove rows when b[i+1] == b[i]
      Iremove <- numeric()
      for (i in 1:(model@n))
      { if (b[i+1] == b[i])
      {  Iremove <- c(Iremove, i)}
      }
      if (length(Iremove) > 0)
      {  b <- b[-Iremove]   
         a <- a[-Iremove]
         b.grad <- b.grad[,-Iremove]
         a.grad <- a.grad[,-Iremove]
      }
      
      # Initialize loop
      nobs <- length(a)-1
      C <- rep(0, nobs+2)
      C[1] <- -1e36
      C[length(C)] <- 1e36
      A1 <- 0
      
      # Loop: build array of indices A1
      for (k in 2:(nobs+1))
      {
        nondom <- 1
        
        if (k == nobs+1)
        {  nondom <- 1
        } else if ( (a[k+1] >= a[k]) && (b[k] == b[k+1]) )
        {  nondom <- 0}
        
        if (nondom == 1)       
        {
          loopdone <- 0
          count <- 0
          while ( loopdone == 0 && count < 1e3 )
          {
            count <- count + 1
            u <- A1[length(A1)] + 1
            C[u+1] <- (a[u]-a[k]) / (b[k] - b[u])
            if ((length(A1) > 1) && (C[u+1] <= C[A1[length(A1)-1]+2]))
            { A1 <- A1[-length(A1)]
            } else
            { A1 <- c(A1, k-1)
              loopdone <- 1
            }
          }
        }
      }
      # Build 'tilde' data
      at <- a[A1+1]
      bt <- b[A1+1]
      ct <- C[c(1, A1+2)]
    } else
    { 
      #--------------- Load data from envir -----------------------------
      Isort <- envir$Isort
      Iremove <- envir$Iremove
      A1 <- envir$A1
      at <- envir$at
      bt <- envir$bt
      ct <- envir$ct
      ######## Sort and reduce A and B ####################
      b.grad <- b.grad[,Isort]  
      a.grad <- -a.grad[,Isort]
      
      if (length(Iremove) > 0)
      {  b.grad <- b.grad[,-Iremove]
         a.grad <- a.grad[,-Iremove]
      }
      
    }
    
    nt <- length(at)
    at.grad <- a.grad[,A1+1]
    bt.grad <- b.grad[,A1+1]
    
    ct.grad <- matrix(0,model@d,nt+1)
    
    for (i in 1:(nt-1))
    {
      ct.grad[,i+1] <- (at.grad[,i]-at.grad[,i+1]) - (at[i]-at[i+1])*(bt.grad[,i+1]-bt.grad[,i])/(bt[i+1]-bt[i])
    }
    
    ######### AGK GRADIENT ##########
    grad.AKG <- matrix(0,model@d,1)     
    for (k in 1:nt)
    {  
      grad.AKG <- grad.AKG + at.grad[,k]*(pnorm(ct[k+1])-pnorm(ct[k])) + bt.grad[,k]*(dnorm(ct[k])-dnorm(ct[k+1])) + 
        ct.grad[,k+1]*dnorm(ct[k+1])*(at[k]+bt[k]*ct[k+1]) - ct.grad[,k]*dnorm(ct[k])*(at[k]+bt[k]*ct[k])
    }
    
    # Compute m_min
    m_min <- min(mk.X)
    m_min <- min(m_min,mk.x)
    I.min <- which.min(c(m_min,mk.x))
    
    m_min.grad <- matrix(0,model@d,1)
    if (I.min == 2)
    {  m_min.grad <- A.grad[,model@n+1]}
    
    grad.AKG <- grad.AKG - (m_min.grad)
  }
  return(grad.AKG)
}
