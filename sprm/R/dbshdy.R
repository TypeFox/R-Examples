dbshdy <-
function(X,y,h,opt="b"){
  
  # % DBSHDY The Jacobian matrix with respect to y in 
  # % univariate SIMPLS regression.
  # % Inputs are:
  # % X, y: datamatrices - supposed to be centered;
  # % h: No. of latent variables;
  # % opt: option, if set to 'b', the Jacobian matrix
  # %      for the regression vector is computed, if 
  # %      set to 'y', the Jacobian for the predicted
  # %      concentrations is computed. 
  # %
  # % Outputs: The jacobian and the regression vector.
  # %
  # % Written by S. Serneels for MiTAC, Uni Antwerp, April 2003.
  
  # Reference: Serneels, S., Van Espen, P.J., Sample specific prediction intervals in SIMPLS, 
  # in: PLS and related methods, M. Vilares, M. Tenenhaus, P. Coelho, V. Esposito Vinzi, A. Morineau (eds.), 
  #     DECISIA, Levallois Perret (France), 2003, pp. 219-233.  
  
  #  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  Unisimpls <- function(X, y, a) {
    n <- nrow(X)
    px <- ncol(X)
    X <- scale(X,center=TRUE,scale=FALSE)
    y <- scale(y,center=TRUE,scale=FALSE)
    if (px > n) {
      dimensions <- 1
      dimension <- px - n
      ressvd <- svd(t(X))
      X <- ressvd$v %*% diag(ressvd$d)
      n <- nrow(X)
      px <- ncol(X)
    }
    else {
      dimensions <- 0
    }
    s <- t(X) %*% y
    U <- matrix(0, nrow = n, ncol = a)
    V <- matrix(0, nrow = px, ncol = a)
    VV <- matrix(0, nrow = px, ncol = a)
    P <- matrix(0, nrow = px, ncol = a)
    R <- matrix(0, nrow = px, ncol = a)
    RR <- matrix(0, nrow = px, ncol = a)
    B <- V
    q <- 1
    for (j in 1:a) {
      r <- s
      RR[, j] <- r
      u <- X %*% r
      # u <- u - U[, 1:max(1, j - 1)] %*% (t(U[, 1:max(1, 
      #                                               j - 1)]) %*% u)
      normu <- drop(sqrt(t(u) %*% u))
      u <- u/normu
      r <- r/normu
      p <- t(X) %*% u
      P [ ,j] <- p
      v <- p - V[, 1:max(1, j - 1)] %*% (t(V[, 1:max(1, j - 1)]) %*% p)
      VV [ ,j] <- v
      v <- v/drop(sqrt(t(v) %*% v))
      s <- s - v %*% (t(v) %*% s)
      q <- drop(t(y)%*%u)
      U[, j] <- u
      R[, j] <- r
      V[, j] <- v
      B[, j] <- R[, 1:j] %*% t(R[, 1:j]) %*% t(X) %*% y
    }
    if (dimensions == 1) {
      B <- ressvd$u %*% B
      P <- ressvd$u %*% P 
      VV <- ressvd$u %*% VV
      R <- ressvd$u %*% R
      RR <- ressvd$u %*% RR
    }
    return(list(scores = U, loadings = P, loadings.deflated = VV, weightings = R, weightings.unscaled = RR))
  }
  
  # --------------------------------------------------------------------------------
  # --------------------------------------------------------------------------------  
  
  n <- nrow(X)
  p <- ncol(X)  
  if(as.integer(h) > as.integer(min(n,p))){
    print(h)
    print(paste("nrow", n))
    print(paste("ncol", p))
    stop('Max # LV cannot exceed the rank of the data matrix')
  }
  X <- scale(X,center=TRUE,scale=TRUE)
  y <- scale(y,center=TRUE,scale=TRUE)  
  S=t(X)%*%X
  s=t(X)%*%y
  res.pls <- Unisimpls(X,y,h)
  
  th <- res.pls$scores[,1]
  b <- matrix(0,nrow=p,ncol=h)
  dsdy <- t(X)
  
  
  for(i in 1:h){
    
    if(i==1){
      
      ah <- res.pls$weightings.unscaled[,1]
      rh <- res.pls$weightings[,1]
      dahdy <- dsdy
      drhdy <- dahdy/as.numeric(sqrt(t(ah)%*%S%*%ah))-(ah%*%t(ah)%*%t(X)%*%X%*%dahdy)/as.numeric(sqrt(t(ah)%*%S%*%ah)^3)
      if(opt=="y"){dthdX <- X%*%drhdy}
      dphdy <- S%*%drhdy
      dvhdy <- dphdy
      if(opt=="y"){
        dyhhdy=((t(y)%*%th%*%diag(n) + kronecker(t(y),th))-2*t(y)%*%th%*%th%*%t(th))%*%dthdy
      }
      else{
        b[,i] <- rh%*%t(rh)%*%s;
        dbhdy <- (kronecker(rh,t(s))+as.numeric(t(rh)%*%s)*diag(p))%*%drhdy + rh%*%t(rh)%*%dsdy
      }
      vh <- res.pls$loadings.deflated[,i]
      
    } else {
      
      dahdy <- (dahdy  - (vh%*%t(vh)/as.numeric(t(vh)%*%vh))%*%dahdy 
                - ((as.numeric(t(ah)%*%vh)*diag(p)+kronecker(t(ah),vh))/as.numeric(t(vh)%*%vh)
                   - (2*as.numeric(t(ah)%*%vh)*vh%*%t(vh))/as.numeric(t(vh)%*%vh)^2)%*%dvhdy)
      
      ah <- res.pls$weightings.unscaled[,i] 
      rh <- res.pls$weightings[,i] 
      drhdy <- dahdy/sqrt(as.numeric(t(ah)%*%S%*%ah)) -(ah%*%t(ah))%*%t(X)%*%X%*%dahdy/as.numeric(sqrt(t(ah)%*%S%*%ah)^3)
      if(opt=="y"){dthdy <- X%*%drhdy + kronecker(t(rh),diag(n))}
      dphdy <- S%*%drhdy
      ph <- res.pls$loadings[,i]
      dvhdy <- (dphdy - (as.numeric(t(ph)%*%vh)*diag(p) + kronecker(t(ph),vh))%*%dvhdy/as.numeric(t(vh)%*%vh)
                + 2*as.numeric(t(ph)%*%vh)*vh%*%t(vh)%*%dvhdy/as.numeric(t(vh)%*%vh)^2)
      dvhdy <- dvhdy - vh%*%t(vh)%*%dphdy/as.numeric(t(vh)%*%vh)
      vh <- res.pls$loadings.deflated[,i]
      if(opt=="y"){
        dyhhdy <- dyhhdy +(as.numeric(t(y)%*%th)*diag(n)+kronecker(t(y),th)-2*as.numeric(t(y)%*%th)%*%th%*%t(th))%*%dthdy
      } else {  
        b[,i]  <- b[,i-1] + rh%*%t(rh)%*%s
        dbhdy <- dbhdy + (kronecker(rh,t(s))+as.numeric(t(rh)%*%s)*diag(p))%*%drhdy + rh%*%t(rh)%*%dsdy
      }
    }
  }
  
  ifelse(opt=="y",output <- dyhhdy,output <- list(dbhdy = dbhdy, b=b))
  return(output)  
}
