# This is file ../spam/R/rmvnorm.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     



# Draw from a multivariate normal:
# (Algorithm 2.3 from Rue and Held, 2005)
rmvnorm.spam <- function(n,mu=rep(0, nrow(Sigma)), Sigma, Rstruct=NULL, ...) {
  # taken from ?chol.spam
  
  if (is(Rstruct,"spam.chol.NgPeyton"))
    cholS <- update.spam.chol.NgPeyton( Rstruct, Sigma, ...)
  else
    cholS <- chol.spam( Sigma,...)
  # cholS is the upper triangular part of the permutated matrix Sigma
  iord <- ordering(cholS, inv=TRUE)

  N <- dim(Sigma)[1]

  R <- as.spam(cholS)
  retval <- as.matrix( ( array(rnorm(n*N),c(n,N)) %*% R)[,iord,drop=F])
     # It is often better to order the sample than the matrix
     # R itself.
  return(sweep(retval, 2, mu, "+"))
}

# Draw from a multivariate normal given a precision matrix:
# (Algorithm 2.4 from Rue and Held, 2005)
rmvnorm.prec <- function(n,mu=rep(0, nrow(Q)), Q, Rstruct=NULL, ...) {

  if (is(Rstruct,"spam.chol.NgPeyton"))
    R <- update.spam.chol.NgPeyton( Rstruct, Q, ...)
  else
    R <- chol(Q,...)
  # R is the upper triangular part of the permutated matrix Sigma

  N <- dim(Q)[1]
  nu <- backsolve(R, array(rnorm(n*N),c(N,n)), k=N)
  return(t(nu+mu))
}


# Draw from the canonical representation of a multivariate normal:
# (Algorithm 2.5 from Rue and Held, 2005)
rmvnorm.canonical <- function(n, b, Q, Rstruct=NULL, ...) {
  N=dim(Q)[1]
  if (is(Rstruct,"spam.chol.NgPeyton"))
    R <- update.spam.chol.NgPeyton( Rstruct, Q, ...)
  else
    R <- chol(Q,...)

  if(is(R,"spam.chol.NgPeyton")){
     mu <- drop(solve.spam( R, b))	
  } else {
     mu <- backsolve( R, forwardsolve( t(R), b), k=N)
  }
  nu <- backsolve(R, array( rnorm(n*N), c(N, n)), k=N)
  return(t(nu+mu))
}



rmvnorm.const <- function (n, mu = rep(0, nrow(Sigma)), Sigma,
                                Rstruct = NULL, 
                                A = array(1, c(1,nrow(Sigma))), a=0, U=NULL,  ...) 
{
    N <- dim(Sigma)[1]
    if (!identical(ncol(A), N)) stop("Incorrect constraint specification")

    if (is(Rstruct, "spam.chol.NgPeyton")) 
        cholS <- update.spam.chol.NgPeyton(Rstruct, Sigma, ...)
    else cholS <- chol.spam(Sigma, ...)
    iord <- ordering(cholS, inv = TRUE)
    N <- dim(Sigma)[1]
    R <- as.spam(cholS)
    x <- sweep( (array(rnorm(n * N), c(n, N)) %*% R)[, iord, drop = F], 2, mu, "+")

    if (is.null(U)){
      V <- backsolve( R, forwardsolve( R, t(A)), k=N)
      W <- A %*% V
      U <- solve(W, t(V))
    }
    correct <- A %*% t(x) - a
    return(x - t( t(U)%*% correct))
}


rmvnorm.prec.const <- function (n, mu = rep(0, nrow(Q)), Q,
                                Rstruct = NULL, 
                                A = array(1, c(1,nrow(Q))), a=0, U=NULL,  ...) 
{
    N = dim(Q)[1]
    if (!identical(ncol(A), N)) stop("Incorrect constraint specification")

    if (is(Rstruct, "spam.chol.NgPeyton")) 
        R <- update.spam.chol.NgPeyton(Rstruct, Q, ...)
    else R <- chol(Q, ...)
    x <- backsolve(R, array(rnorm(n * N), c(N, n)), k=N) + mu

       
    if (is.null(U)){
        tV <- t( backsolve( R, forwardsolve( R, t(A)), k=N))
        W <- tcrossprod(A, tV)
        U <- solve(W, tV)
    }
    correct <- A %*% x - a
    return(t(x- t(U) %*% correct))
}



rmvnorm.canonical.const <- function (n, b, Q, Rstruct = NULL, 
                                     A = array(1, c(1,nrow(Q))), a=0, U=NULL, ...) 
{
    N = dim(Q)[1]
    if (!identical(ncol(A), N)) stop("Incorrect constraint specification")

    if (is(Rstruct, "spam.chol.NgPeyton")) 
        R <- update.spam.chol.NgPeyton(Rstruct, Q, ...)
    else R <- chol(Q, ...)
    if (is(R, "spam.chol.NgPeyton")) {
        mu <- drop(solve.spam(R, b))
    }
    else {
        mu <- backsolve(R, forwardsolve(t(R), b))
    }
    x <- backsolve(R, array(rnorm(n * N), c(N, n)), k=N) + mu

   
    if (is.null(U)){
        tV <- t( backsolve( R, forwardsolve( R, t(A)), k=N))
        W <- tcrossprod(A, tV)
        U <- solve(W, tV)
    }
    correct <- A %*% x - a
    return(t(x- t(U) %*% correct))
 
    
}
