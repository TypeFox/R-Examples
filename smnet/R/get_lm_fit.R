
get_lm_fit<-function(P.list, X.spam, XTX.spam, X.list, response, 
                     n.linear = n.linear, sm.names, lin.names, 
                     lin.means){
  
  n         <- length(response)
  X.dim     <- lapply(X.list, ncol)
  n.terms   <- length(X.list)
  np        <- sum(unlist(X.dim))
  inds      <- unlist(X.dim)
  cum.inds  <- cumsum(inds)
  Xy        <- t(X.spam)%*%response
  
  # fit the model
  U          <- chol.spam(XTX.spam)
  beta_hat   <- backsolve.spam(U, forwardsolve.spam(U, Xy))   
  vec        <- forwardsolve.spam(U, t(X.spam))
  dof.params <- rowSums(vec^2)
  dof.list   <- vector("list", length = n.terms)
  dof.list[[1]]<-dof.params[1]
  if(n.terms > 1)  for(i in 2:n.terms) dof.list[[i]] <- dof.params[(cum.inds[i-1]+1):(cum.inds[i])]
  ED       <- sum(dof.params)
  fit      <- X.spam %*% beta_hat
  sigma.sq <- sum((response - fit)^2)/(n - ED - 1) # model variance
  # adjust the model intercept for covariate mean centering
  if(n.linear > 0) beta_hat[1] <- beta_hat[1] - sum(lin.means*beta_hat[2:(1 + n.linear)])
  list(U = U, beta_hat=beta_hat, ED = ED, fit = fit, sigma.sq = sigma.sq, 
       lin.means = lin.means, n.linear = n.linear)
}


