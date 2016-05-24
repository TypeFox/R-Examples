homoStat <- function(y, v) {
  if (is.vector(y)) no.y <- 1 else no.y <- ncol(y)  
  if (is.vector(v)) no.v <- 1 else no.v <- ncol(v)
  if ( no.v != no.y*(no.y+1)/2 )
    stop(paste("The expected no. of columns in v is ", no.y*(no.y+1)/2,
               " while the observed no. of columns in v is ", no.v, ".", sep=""))
    
  if (no.y==1) {
    miss.index <- is.na(y)
    y <- y[!miss.index]
    v <- v[!miss.index]
    w <- 1/v
    beta <- sum(y*w)/sum(w)
    Q <- sum( w*(y-beta)^2 )
    Q.df <- length(y)-1
    pval <- 1-pchisq(Q, df=Q.df)
  } else {
    Y <- matrix( c(t(y)), ncol=1 )
    miss.index <- is.na(Y)
    Y <- matrix( Y[!miss.index], ncol=1 )
    X <- matrix( rep(Diag(no.y), nrow(y)), ncol=no.y, byrow=TRUE )

    X <- X[!miss.index, , drop=FALSE]
    V <- matrix2bdiag(v)
    V <- V[!miss.index, !miss.index, drop=FALSE]
    
    V_inv <- solve(V)
    Q <- t(Y) %*% ( V_inv - V_inv %*% X %*% solve(t(X)
              %*% V_inv %*% X) %*% t(X) %*% V_inv ) %*% Y
    Q.df <- nrow(X)-ncol(X)
	pval <- 1-pchisq(Q, df=Q.df)
  }
    list(Q=Q, Q.df=Q.df, pval=pval)
}


