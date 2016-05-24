EMMA <- function (y, X=NULL, Z=NULL, K=NULL, REML=TRUE, silent=FALSE) { 
  if(!silent){
    count <- 0
    tot <- 1
    pb <- txtProgressBar(style = 3)
    setTxtProgressBar(pb, 0)
  }
  y.or <- y
  ## to be used for fitted values
  if(is.null(X)){ # only random effects present
    x.or <- matrix(rep(1, length(y)), ncol=1)
  }else{x.or <- X}
  if(is.null(Z)){ # only random effects present
    z.or <- diag(length(y))
  }else{z.or <- Z}
  ##
  good <- which(!is.na(y))
  y <- y[good]
  Z <- Z[good,]
  if(is.null(X)){ # only random effects present
    X <- matrix(rep(1, length(y)), ncol=1)
  }else{
    X <- as.matrix(X[good,])
  }
  
  if(!is.null(X)){ # both present, extract xm from X, check double list
    if(is.list(X)){
      if(is.list(X[[1]])){
        X=X[[1]][[1]]
      }else{
        X=X[[1]]
      }
    }else{
      X=X 
    }
  }
  ##############
  if(is.null(Z)){ # only random effects present
    Z <- diag(length(y))
  }
  if(!is.null(Z)){ # both present, extract xm from X, check double list
    if(is.list(Z)){
      if(is.list(Z[[1]])){
        Z=Z[[1]][[1]]
      }else{
        Z=Z[[1]]
      }
    }else{
      Z=Z 
    }
  }
  ##############
  if(is.null(K)){ # only random effects present
    K <- diag(length(y))
  }
  if(!is.null(K)){ # both present, extract xm from X, check double list
    if(is.list(K)){
      if(is.list(K[[1]])){
        K=K[[1]][[1]]
      }else{
        K=K[[1]]
      }
    }else{
      K=K 
    }
  }
  #diag(K) <- diag(K) + 1e-06
  ##-------------------- S = I - X(X'X)-X' ---------------------##
  ## usually is used to estimate residuals as (sometimes called projection matrix): 
  ## y - y.hat = y-Xb = y - X(X'X)X'y = [I - X(X'X)-X']y = e = residuals
  ## I - X(X'X)-X' is symmetric and idempotent, where X(X'X)-X' 
  ## is the Hat matrix to project y on X or residuals on y
  n = length(y) # number of observations
  p <- dim(X)[2] # number of columns
  II <- as(diag(n), "sparseMatrix")
  S <- diag(n) -  X %*% crossprod(  solve(crossprod(X))  , t(X)  )
  S <- as(S, "sparseMatrix")
  K <- as(K, "sparseMatrix")
  Z <- as(Z, "sparseMatrix")
  ##-------------------- H = ZKZ' + delta.prov*I ---------------------##
  de <- log(n) #initial value of delta [Var(e)/Var(u)]
  H <- (  Z %*% (K%*%t(Z)) ) + ( de * II )
  ## ------------------- eigen(SHS) --> UV ---------------------##
  eigSHS <- eigen(S %*% H %*% S, symmetric = TRUE)
  lambda <- eigSHS$values[1:(n-p)] - de # # eigen values of SHS = [I - X(X'X)-X'] [ZKZ' + I*theta] [I - X(X'X)-X']
  U <- eigSHS$vectors[,(1:(n-p))]    # eigen vectors of SHS = [I - X(X'X)-X'] [ZKZ' + I*theta] [I - X(X'X)-X']
  # if ML
  Hb.system <- eigen(H, symmetric = TRUE)
  phi <- Hb.system$values - de
  # rotation matrix that doesn't change the direction of the data 
  # eigenvalues scale the vectors of the rotation matrix
  # is the same than having (ZKZ+Ide)-
  
  ##-------------------- estimate parameter eta and delta ---------------------## 
  # part of "y" variability attributed to random effects 
  eta <- crossprod(U,y) # t(U) %*% y   Uy ; where the eigenvectors of random effects project weighted by A matrix 
  # y / (ZKZ+Ide)
  
  ##-------------------- likelihood function using eta, lambda and delta ---------------------##
  # eta = y / (ZKZ+Ide)
  # lambda= eigenvalues((ZKZ+Ide)-)
  # delta = Var(e)/Var(u)
  if(REML == TRUE){
    minimfunc <- function(delta) {  (n - p) * log(sum(eta^2/(lambda + delta))) + sum(log(lambda + delta))   }
  }else{
    minimfunc <- function(delta) {n * log(sum(eta^2/(lambda + delta))) + sum(log(phi + delta)) }
  }
  ##-------------------- REML estimator of delta=Var(e)/Var(u)---------------------##
  REML <- optimize(minimfunc, lower = 9^(-9), upper = 9^9, tol = 1e-06)
  delta <- REML$minimum
  ##-------------------- H- = (H = ZKZ' + delta*I)- ---------------------## 
  # remember V- = 1/var(g) * H-  which solves the random effect model, here there's an extra ridge parameter
  hhhh <- Z%*%K%*%t(Z) +  delta*diag(length(y)) 
  H.hat.inv <- solve(as(hhhh ,"sparseMatrix"))
  ##-------------------- FIXED EFFECTS Beta --> X'H-X = X'H-y  -->   B = (X'H-X)-X'Hy 
  # regular linear model using the weighted matrix H- which takes into account random effects and correlation matrices
  aaa <- t(X)%*%H.hat.inv%*%X
  bbb <- t(X)%*%H.hat.inv%*%y
  beta <- solve(aaa, bbb)
  error <- y - (X%*%beta) # residuals e = y-XB
  #He <- H.hat.inv %*% error  # e / (ZKZ+Ide) , part of error due to random effects
  ##-------------------- Variance components Var(u) and Var(e) ---------------------##
  # SSu / (n-p)
  sigma2.u <- sum ( (eta^2) / (lambda+delta) )   /   (n - p)
  sigma2.e <- delta * sigma2.u  # Var(e)/Var(u) * Var(u) = Var(e) !!!!
  ##-------------------- RANDOM EFFECTS u = (ZK)'H- ---------------------##
  u <- t( Z %*% K ) %*% (H.hat.inv %*% error)  # ZK H- (y-XB) # gianola uses letter V- instead of H-
  df <- n - p
  pi <- 3.14159
  ll <- -0.5 * (  REML$objective + df + ( df * log((2*pi)/df) )  )
  ## NOW we will focus in estimating variances and PEV for BLUP's called u
  ##-------------------- V inverse is used for estimating P and Var(u)  ---------------------##
  # [1/Var(u)] * H-   is the variance for random effects
  Vi <- (1/sigma2.u) * H.hat.inv
  ##-------------------- Vi  -  ViX  solve(X'VX = X'V) ---------------------##
  P <- Vi - (    (Vi%*%X) %*%  solve(crossprod(X, Vi %*% X), crossprod(X, Vi))    )#; dim(P)
  ##-------------------- Variance u, B, ---------------------##
  # variance for the variance or 2nd moment sigma2^2 = sigma4
  Var.u <- (sigma2.u^2) *  ( crossprod(Z%*%K, P)  %*%  (Z%*%K)   ) # sigma^4 ZKP ZK
  PEV.u <- sigma2.u * K - Var.u  # standard errors (SE) for each individual
  var.beta <- solve(crossprod(X, Vi %*% X)) # variance of fixed effects
  
  out <- matrix(c(sigma2.u,sigma2.e))
  rownames(out) <- c("V(u)","V(e)")
  
  AIC = (-2 * ll ) + ( 2 * dim(X)[2]) # k=2, npar=2
  BIC = (-2 * ll ) + ( log(length(y)) * dim(X)[2])
  
  fitted.y <- x.or %*% beta
  fitted.y <- fitted.y + (z.or %*% u)
  fitted.u <- (z.or %*% u)
  
  fitted.y.good <- fitted.y[good]
  residuals2 <- y.or[good] - fitted.y[good] # conditional residuals
  
  rownames(u) <- colnames(Z)
  rownames(beta) <- colnames(X)
  
  if(!silent){
    count <- count + 1
    setTxtProgressBar(pb, (tot/tot))
  }
  return(list(var.comp=out, beta.hat = beta, 
              u.hat = u, Var.u.hat = (Var.u), 
              Var.beta.hat = var.beta, PEV.u.hat = (PEV.u), 
              LL = ll, AIC=AIC, BIC=BIC, V.inv=H.hat.inv, X=X, Z=Z, K=K, 
              fitted.y=fitted.y, fitted.u=fitted.u, residuals=error, 
              cond.residuals=residuals2, fitted.y.good=fitted.y.good))
  # diag(Var(u)) is variance for random effects, each marker or plant
  # diag(Var(beta))
}

