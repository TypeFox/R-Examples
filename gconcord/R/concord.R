#' CONvex CORrelation selection methoD
#'
#' Estimates a sparse inverse covariance matrix from a convex 
#' pseudo-likelihood function with lasso L1 penalty
#'
#' Implements the CONCORD method by Khare, Oh and Rajaratnam (2013)
#' http://arxiv.org/abs/1307.5381
#'
#' @param data Data matrix with n observations (rows) and
#'        p variables (columns)
#' @param lambda Penalty parameter
#' @param tol Convergence threshold
#' @param maxit Maximum number of iterations before termination
#' @param save.iterates Returns iterates if TRUE
#' @param ... ignored
#'
#' @useDynLib gconcord
#'
#' @examples
#' library(mvtnorm)
#'
#' ## True omega
#' omega <- matrix(0,3,3)
#' omega[1,2] <- omega[2,3] <- 2.1
#' omega <- t(omega) + omega
#' diag(omega) <- 3
#' 
#' sigma <- solve(omega)
#' 
#' ## Generate data
#' set.seed(60)
#' data <- rmvnorm(100, rep(0,3), sigma)
#' 
#' ## Solve
#' concord(data,2)
#'
#' @export
concord <- function(data, lambda, tol=1e-5, maxit=100, save.iterates=FALSE, ...){

  n <- nrow(data)
  p <- ncol(data)

  omega <- diag(p)  

  r <- ifelse(save.iterates, maxit, 1)
  iterates <- matrix(0, r, p*(p+1)/2)
  iterates[1] <- ifelse(save.iterates, 1, 0) ## flag for concordC to save iterates

  olambda <- lambda
  lambda <- 2*lambda
  if (length(lambda)==p^2) {
    lambda <- lambda[lower.tri(lambda,diag=TRUE)]
  } else if (length(lambda)==1) {
    lambda <- rep(lambda, p*(p+1)/2)
  } else {
    stop("wrong number of elements in 'lambda'")
  }

  if (n<=p) { 
      S <- diag(colSums(data*data)/n)
      out <- .C("concordC2",
                n=as.integer(n),
                p=as.integer(p),
                Y=as.double(data),
                S=as.double(S),
                lambda=as.double(lambda),
                omega=as.double(omega),
                tol=as.double(tol),
                maxit=as.integer(maxit),
                iterates=as.double(iterates),
                PACKAGE='gconcord')
  } else {
      S <- (t(data)%*%data)/n
      out <- .C("concordC",
                n=as.integer(n),
                p=as.integer(p),
                S=as.double(S),
                lambda=as.double(lambda),
                omega=as.double(omega),
                tol=as.double(tol),
                maxit=as.integer(maxit),
                iterates=as.double(iterates),
                PACKAGE='gconcord')
  }
  out$omega <- matrix(out$omega,p,p)

  output <- structure(out$omega, lambda=olambda, iter=out$maxit,
                      solver='concord', tol=tol)
                      
  if (save.iterates) {
    out$iterates <- matrix(out$iterates,r,p*(p+1)/2, byrow=TRUE)

    if (out$maxit<maxit){
      out$iterates <- out$iterates[1:out$maxit,]
    }
    
    attr(output,'iterates') <- out$iterates
  }
  
  output
}

#' Symmetric Lasso (symlasso)
#'
#' Estimates a sparse inverse covariance matrix from a
#' pseudo-likelihood function formulation with L1 penalty
#' on inverse covariance elements.
#' 
#' Implements the Symmetric Lasso method by Friedman, 
#' Hastie and Tibshirani (2010) 
#' http://statweb.stanford.edu/~tibs/ftp/ggraph.pdf
#'
#' @param data Data matrix with n observations (rows) and
#'        p variables (columns)
#' @param lambda Penalty parameter
#' @param tol Convergence threshold
#' @param maxit Maximum number of iterations before termination
#' @param save.iterates Returns iterates if TRUE
#' @param ... ignored
#'
#' @examples
#' library(mvtnorm)
#' 
#' ## True omega
#' omega <- matrix(0,3,3)
#' omega[1,2] <- omega[2,3] <- 2.1
#' omega <- t(omega) + omega
#' diag(omega) <- 3
#' 
#' sigma <- solve(omega)
#' 
#' ## Generate data
#' set.seed(60)
#' data <- rmvnorm(100, rep(0,3), sigma)
#' 
#' ## Solve
#' symlasso(data,2.1)
#'
#' @export
symlasso <- function(data, lambda, tol=1e-5, maxit=100, save.iterates=FALSE, ...){

  p <- ncol(data)
  n <- nrow(data)
  
  Cmat <- (t(data)%*%data)/n
  omega <- matrix(0,p,p)
  
  r <- ifelse(save.iterates, maxit+1, 1)
  iterates <- matrix(0, r, p*(p+1)/2)
  iterates[1] <- ifelse(save.iterates, 1, 0) ## flag for symlassoC to save iterates
  
  out <- .C("symlassoC",
            n=as.integer(n),
            p=as.integer(p),
            Cmat=as.double(Cmat),
            lambda=as.double(lambda),
            omega=as.double(omega),
            tol=as.double(tol),
            maxit=as.integer(maxit),
            iterates=as.double(iterates),
            PACKAGE='gconcord')

  out$omega <- matrix(out$omega,nrow(Cmat),ncol(Cmat))

  output <- structure(out$omega, lambda=out$lambda, iter=out$maxit,
                      solver='symlasso', tol=tol)

  if (save.iterates) {
    out$iterates <- matrix(out$iterates,r,p*(p+1)/2, byrow=TRUE)

    if (out$maxit<maxit){
      out$iterates <- out$iterates[1:out$maxit,]
    }
    
    attr(output,'iterates') <- out$iterates
  }
  
  output
}
