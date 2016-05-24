#@title Function that computes a submatrix of metric M^(-1/2)
#@export 
#@param x the column vector of a regressor or factor
#' @importFrom stats model.matrix
metricBloc <- function(x){
  z <- model.matrix(~x)[,-1,drop=FALSE]
  z <- scale(z,center=TRUE,scale=FALSE)
  z <- (t(z)%*%z)/nrow(z)
  return(z)  
}


#@title Function that computes M^(-1/2)
#@export 
#@param data data.frame containing all covariates
metric <- function(data){
  z <- sapply(data,metricBloc)
  z <- bdiag(z)
  if(ncol(z)==1){
    z <- diag(as.vector(as.matrix(z)))
  }
  z <- svd(z)
  z <- z$u%*%diag(1/sqrt(z$d))%*%t(z$v)
  return(z)
}
