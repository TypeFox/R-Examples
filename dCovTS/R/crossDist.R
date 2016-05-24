crossDist <- function(x,lags){
 n <- NROW(x)
 p <- NCOL(x)
 Y <- x[1:(n-lags),]
 X <- x[(1+lags):n,]
 A <- lapply(1:p,function(j) as.matrix(dist(X[,j],diag=TRUE,upper=TRUE)))
 B <- lapply(1:p,function(j) as.matrix(dist(Y[,j],diag=TRUE,upper=TRUE)))
 return(list(A=A,B=B))
}


