#### Sandwich Estimator

### vov.sandwich.glm caluclates the sandwich esitmator of the
### standard errors obtained by 'glm'.
###
### Author: Zhushan "Mandy" Li
### Last Updated: July 29, 2006


### version 3
### This version is kind of compromised comparing to the original
### version(at the end of the file, commented out by if(0){...} ). The
### function limits that the strata variable must strictly follow the
### format of 'cid' given in plStackData.
vcov.sandwich.glm <- function(object, adjust=FALSE, strata=NULL, ...){
  X <- model.matrix(object)

  wt <- weights(object, "prior")
  r <- residuals(object, "response")

  v <- r * sqrt(wt)
  Xv <- v * X  # row mulitply vector r*sqrt(wt) 

  if(is.null(strata)){
    J <- crossprod(Xv)
  } else {
    ## it is required strata be in the strict format of (1, 2, 3, ...,
    ## n, 1,2,3,...,n, 1,2,3,...,n,.....), otherwise the following
    ## code does not work.
    nlen <- length(strata)
    nstrat <- strata[nlen]
    nrep <- nlen/nstrat
    jjstrata <- rep(1:nstrat, nrep)
    if(any(jjstrata != strata)){ stop('strata is not as expected') }

    ncolXv <- ncol(Xv)
    dim(Xv) <- c(nstrat, nrep, ncolXv)
    Xvv <- apply(Xv, c(1,3), sum) 
    
    J <- crossprod(Xvv)
  }
  
  covb <- vcov(object)
  covb.sandwich <- covb %*% J %*% covb

  if(adjust){
    n <- sum(wt)
    p <- object$rank
    covb.sandwich <- n/(n-p) * covb.sandwich
  }
  
  return (covb.sandwich)
}


if(0){ ## version 2
#### Trying to address the not-enough-memory problem, I used the
#### sparse matrix in package 'Matrix', however the for loop in the
#### function is too slow for large data set
vcov.sandwich.glm <- function(object, adjust=FALSE, strata=NULL){
  X <- Matrix(model.matrix(object))

  wt <- weights(object, "prior")
  r <- residuals(object, "response")

  if(is.null(strata)){
#    v <- r * sqrt(wt)
#    Xv <- sweep(X, 1, v, '*')
#    J <- crossprod(Xv)
    V <-  Diagonal(x=r^2 * wt)
  } else {
#    strata.mtx <- outer(strata, strata, '==')
#    V <- outer(r,r,'*') * (diag(wt) %*% strata.mtx)    
    n <- nrow(X)
#    strata.mtx <- as(Matrix(0, n, n), "sparseMatrix")
    V<- as(Matrix(0, n, n), "sparseMatrix")
    rw <- r*sqrt(wt)
    for(s in unique(strata)){
      spos <- strata==s
#      strata.mtx[spos, spos] <- 1
      V[spos, spos] <- as.vector(outer(rw[spos], rw[spos], '*'))
    }
  }
  
  J <- t(X) %*% (V %*% X)

  covb <- vcov(object)
  covb.sandwich <- covb %*% J %*% covb

  if(adjust){
    n <- sum(wt)
    p <- object$rank
    covb.sandwich <- n/(n-p) * covb.sandwich
  }
  
  return (covb.sandwich)
}

}

if(0){  ## original version
#### The following code works well for moderate size data. However it
#### does not work for large data set because the matrix V is too
#### large for R to allocate enough memory for it. On the other hand,
#### V is a (block) diagonal matrix so that it usually is sparse. A
#### better and effcient algorithm can be used.
vcov.sandwich.glm <- function(object, adjust=FALSE, strata=NULL){
#  X <- model.matrix(formula(object), data=object$data)
  X <- model.matrix(object)
#  W <- diag(object$weights)
#  Info <-  t(X) %*% W %*% X
  
  wt <- weights(object, "prior")
  r <- residuals(object, "response")

##   if(adjust){
##     ## using leverage to adjust, get an unbiased sandwich estimator.
##     ## Ref: Carroll, Wang, Simpson, STromberg & Ruppert, 1998, The Sandwich
##     ## (Robust Covariance Matrix) Esitmator. 
##     hat.matrix <- sqrt(W) %*% X %*% solve(t(X)%*%W%*%X) %*% t(X) %*% sqrt(W)
##     wt <- wt/(1-diag(hat.matrix))
##   }
  
  if(is.null(strata)){
    V <- diag(r^2 * wt)
  } else {
    strata.mtx <- outer(strata, strata, '==')
    V <- outer(r,r,'*') * (diag(wt) %*% strata.mtx)    
  }
  
  J <- t(X) %*% V %*% X

  covb <- vcov(object)
  covb.sandwich <- covb %*% J %*% covb

  
##   if(adjust){
#### Ref: Hardin & Hilbe, 2003, Generalized Estimating Equations,pp31
##     if(is.null(strata)){
##       n <- nrow(V)
##       p <- object$rank
##     covb.sandwich <- n/(n-p) * covb.sandwich    
##     } else {
##       n <- length(unique(strata))
##     covb.sandwich <- n/(n-1) * covb.sandwich    
##     }
##   }
  if(adjust){
    n <- sum(wt)
    p <- object$rank
    covb.sandwich <- n/(n-p) * covb.sandwich
  }
  
  return (covb.sandwich)
}
  
}
