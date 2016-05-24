edge.prob <-
function(W, log = TRUE, account.prior = FALSE, q0 = 0.5){
  p <- nrow(W)
  if (log){
    M <- exp(W - mean(W[upper.tri(W)]))
    diag(M) <- 0
  } else {
    M <- W
  }
  Delta <- -M+diag(apply(M,1,sum))
  
  prob <- matrix(0,p,p)
  
  Q <- inv_RcppEigen(Delta[-1,-1])
  P <- sapply(1:(p-1), function(x) -2*Q[, x] + Q[x,x])
  P <- t(sapply(1:(p-1), function(x) P[x, ] + Q[x,x]))
  
  prob[-1,-1] <- P
  prob[1,-1] <- diag(Q)
  prob[-1,1] <- diag(Q)
  
  prob <- prob*M

  if (account.prior){return(account.for.prior(prob,q0))} else {return(prob)}
}
