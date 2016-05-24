get_max_df <- function(P.flat, rhoInput, ridgeM, XTX.spam, cholFactor, info, Xy, X.spam, X.list, max.df){
  X.dim     <- lapply(X.list, ncol)
  n.terms   <- length(X.list)
  np        <- sum(unlist(X.dim))
  inds      <- unlist(X.dim)
  cum.inds  <- cumsum(inds) 
  P         <- P.flat
  get_df <- function(rho, rhoInput, ridgeM, X.spam, cholFactor, XTX.spam, P, Xy, n.terms, cum.inds){
    rhoInput[length(rhoInput) - 1] <- rho
    for(j in 1:length(rhoInput)) P[[j]]<-P[[j]]*exp(rhoInput[j])
    info     <- XTX.spam + Reduce("+", P) + ridgeM
    U        <- update.spam.chol.NgPeyton(cholFactor, info)
    beta_hat <- backsolve.spam(U, forwardsolve.spam(U, Xy))  
    left1    <- forwardsolve.spam(U, t(X.spam))
    trH      <- rowSums(left1*left1)
    abs(sum(trH[(cum.inds[n.terms-1]+1):(cum.inds[n.terms])]) - max.df)
  }
  opt_df <- optimize(get_df, rhoInput = rhoInput, interval = c(-20, 20), P = P, ridgeM = ridgeM, XTX.spam = XTX.spam, 
                     cholFactor = cholFactor, Xy = Xy, n.terms = n.terms, cum.inds = cum.inds, X.spam = X.spam)
  opt_df$minimum
}
