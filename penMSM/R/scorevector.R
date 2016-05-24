scorevector <- function(beta, X, risksetlist, event){
  n <- length(event)
  f <- as.numeric(X %*% beta)
  ef <- exp(f)
  p <- length(beta)
  riskmatrix <- matrix(nrow = n, ncol = p, 0)
  for (i in 1:n){
    riskset <- risksetlist[[i]]
    ef.riskset <- ef[riskset]
    currentrisk <- sum(ef.riskset)
    X.i <- X[riskset, ]/currentrisk
    riskmatrix[i, ] <- t(ef.riskset) %*% X.i
    }
  scorevector <- as.numeric(event %*% (X - riskmatrix))
  return(scorevector)}