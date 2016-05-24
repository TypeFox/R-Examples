##
## Marginal Contribution to Risk of a Portfolio
##
mrc <- function(weights, Sigma, percentage = TRUE){
  Sigma <- as.matrix(Sigma)
  if(!isSymmetric(Sigma)){
    stop("Object provided for 'Sigma' is not a symmetric matrix.\n")
  }
  w <- as.vector(weights)
  if(length(w) != ncol(Sigma)){
    stop("Length of 'weights' vector differs from row/column dimension of 'Sigma'.\n")
  }
  sigma <- c(sqrt(t(w) %*% Sigma %*% w))
  sw <- Sigma %*% w
  dw <- c(w * sw/sigma)
  ifelse(percentage, res <- dw / sum(dw) * 100, res <- dw)
  return(res)
}
