getModelComponents <-
function(m, analytic) { 
  # A function that calculates all components needed to calculate the bias
  # correction as in Greven & Kneib (2010)
  #
  # Args: 
  #   m     = Object of class lmerMod. Obtained by lmer()
  #   analytic = FALSE if the numeric hessian of the (restricted) marginal log-
  #              likelihood from the lmer optimization procedure should be used.
  #              Otherwise (default) TRUE, i.e. use a analytical version that 
  #              has to be computed.
  #
  # Returns:
  #   model = List of components needed to calculate the bias correction
  #   
  model         <- list()
  model$df      <- NULL
  X             <- getME(m, "X")
  n             <- nrow(X)
  Z             <- getME(m, "Z")
  theta         <- getME(m, "theta")
  Lambda        <- getME(m, "Lambda")
  Lambdat       <- getME(m, "Lambdat")
  model$Wlist   <- list()
  model$eWelist <- list()
  V0inv         <- diag(rep(1, n)) - crossprod(solve(getME(m, "L"), system = "L") %*% 
                   solve(getME(m, "L"), Lambdat, system = "P") %*% t(Z))
# P             <- diag(rep(1, n)) - X %*%  chol2inv(getME(m, "RX")) %*% crossprod(X, V0inv)
  
  ## pre calculate matrices for faster computation
# A   <- crossprod(P, V0inv)
  A   <- V0inv - crossprod(crossprod(X %*% solve(getME(m, "RX")), V0inv))
  y   <- getME(m, "y") 
  e   <- y - getME(m, "mu")
  
  ## prepare list of derivative matrices W_j
  ind <- getME(m, "Lind")
  len <- rep(0, length(Lambda@x))
  
  for(s in 1:length(theta)) {
#  model$Wlist <- lapply(theta, function(s){
    LambdaS                    <- Lambda
    LambdaSt                   <- Lambdat
    LambdaS@x                  <- LambdaSt@x                    <- len
    LambdaS@x[which(ind == s)] <- LambdaSt@x[which(ind == s)]   <- 1
    diagonal                   <- diag(LambdaS)
    diag(LambdaS)              <- diag(LambdaSt)                <- 0
    Ds                         <- LambdaS + LambdaSt
    diag(Ds)                   <- diagonal
    model$Wlist[[s]]           <- Z %*% Ds %*% t(Z)
    model$eWelist[[s]]         <- as.numeric(e %*% model$Wlist[[s]] %*% e)
#   model$Wlist[[s]]  <- model$Wlist[[s]]/norm(model$Wlist[[s]], type = "F")
  }
  
  ## Write everything into a return list
  model$X       <- X        
  model$n       <- n        
  model$theta   <- theta    
  model$Z       <- Z      
  model$Lambda  <- Lambda 
  model$Lambdat <- Lambdat
  model$V0inv   <- V0inv
  model$A       <- A
  if(analytic) {
    model$B     <- matrix(0, length(theta), length(theta)) 
  } else {
    model$B     <- m@optinfo$derivs$Hessian  
  }
  model$C       <- matrix(0, length(theta), n)   
  model$y       <- y   
  model$e       <- e  
  model$tye     <- as.numeric(crossprod(y, e))
  model$isREML  <- isREML(m)
  
  return(model)
}
