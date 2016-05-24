####################################################################
## FILE CONTAINING DERIVATIVE COMPUTATIONS FOR THE LOGLIKELIHOODS ##
####################################################################
##Functions in this file:
## genGradient            EX:ok
## genHessian             EX:with genGradient
## loglikeSTGrad          EX:ok-dont run
## loglikeSTHessian       EX:with loglikeSTGrad
## loglikeSTnaiveGrad     EX:with loglikeSTGrad
## loglikeSTnaiveHessian  EX:with loglikeSTGrad

##' Computes finite difference gradient and/or hessian. \code{genGradient}
##' function does forward, backward or central differences, the
##' \code{genHessian} function uses only central differences.
##'
##' @title Compute Finite Difference Gradient and Hessians.
##' 
##' @param x Point at which to compute the gradient or hessian.
##' @param func function that takes only \code{x} as an input argument. Use \cr
##'   \code{function(x){my.func(x,other.input)}} to create a temporary
##'   function, see the example.
##' @param h Step length for the finite difference.
##' @param diff.type Type of finite difference, \code{diff.type>0} gives forward
##'   differences, \code{diff.type=0} gives central differences, and
##'   \code{diff.type<0} gives backward differences.
##' 
##' @return gradient vector or Hessian matrix.
##' 
##' @example Rd_examples/Ex_genGradient.R
##' 
##' @author Johan Lindström
##' 
##' @family numerical derivatives
##' @export
genGradient <- function(x, func, h=1e-3, diff.type=0){
  diff.type <- sign(diff.type) ##ensure diff.type is [-1,0,1]
  if(diff.type==0){ ##central
    f.p <- double(length(x))
    f.m <- double(length(x))
    for(i in 1:length(x)){
      dx <- x
      dx[i] <- x[i]+h/2
      f.p[i] <- func(dx)
      dx[i] <- x[i]-h/2
      f.m[i] <- func(dx)
    }
    df <- (f.p-f.m)/h
  }else{ ##forward or backward
    f <- func(x)
    df <- double(length(x))
    for(i in 1:length(x)){
      dx <- x
      dx[i] <- dx[i]+diff.type*h
      df[i] <- func(dx)
    }
    df <- diff.type*(df-f)/h
  }
  return(df)
}##function genGradient

##' @rdname genGradient
##' @export
genHessian <- function(x, func, h=1e-3){
  ##uses only central differences
  N <- length(x)
  ##allocate matrix for the hessian
  hessian <- matrix(NA,N,N)
  ##start with the diagonal elements
  f <- double(5)
  for(i in 1:N){
    tmp <- x
    for(j in 1:5){
      tmp[i] <- x[i]+(j-3)*h
      f[j] <- func(tmp)
    }
    hessian[i,i] <- (-f[1]+16*f[2]-30*f[3]+16*f[4]-f[5])/(12*h^2)
  }
  ##off diagonal elements
  f <- double(4)
  Ind <- cbind(c(1,1,-1,-1),c(1,-1,1,-1))
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      tmp <- x
      for(k in 1:4){
        tmp[i] <- x[i] + Ind[k,1]*h
        tmp[j] <- x[j] + Ind[k,2]*h
        f[k] <- func(tmp)
      }
      hessian[i,j] <- (f[1]-f[2]-f[3]+f[4])/(4*h^2)
      ##use symmetri for the other half of the matrix
      hessian[j,i] <- hessian[i,j]
    }
  }
  return(hessian)
}##function genHessian

##' Computes finite difference gradients and hessians for the log-likelihood
##' functions \code{\link{loglikeST}} and \code{\link{loglikeSTnaive}}.
##' \cr
##' Uses \code{\link{genGradient}} and \code{\link{genHessian}} to compute
##' finite difference derivatives of the log-likelihood function in
##' \code{\link{loglikeST}} and \code{\link{loglikeSTnaive}}.
##' 
##' @title Compute Gradient and Hessian for the Log-likelihood
##' 
##' @param x Point at which to compute the gradient or hessian, see
##'   \code{\link{loglikeST}}.
##' @param STmodel \code{STmodel} object with the model for which to compute
##'   derivatives of the log-likelihood.
##' @param type A single character indicating the type of log-likelihood to
##'   compute. Valid options are "f", "p", and "r", for \emph{full},
##'   \emph{profile} or \emph{restricted maximum likelihood} (REML).
##' @param x.fixed Parameters to keep fixed, see \code{\link{loglikeST}}.
##' @param h,diff.type Step length and type of finite difference to use when
##'   computing gradients, see \code{\link{genGradient}}.
##' 
##' @return Returns the gradient or Hessian for the \code{\link{loglikeST}}
##'   and \code{\link{loglikeSTnaive}} functions.
##' 
##' @section Warning: \code{loglikeSTnaiveGrad} and
##'   \code{loglikeSTnaiveHhessian} may take \strong{very} long time to run,
##'   use with \strong{extreme caution}.
##' 
##' @example Rd_examples/Ex_loglikeSTGrad.R
##' 
##' @author Johan Lindström
##' 
##' @family likelihood functions
##' @family numerical derivatives
##' @export
loglikeSTGrad <- function(x, STmodel, type="p", x.fixed=NULL,
                         h=1e-3, diff.type=0){
  func <- function(x0){ loglikeST(x0, STmodel, type, x.fixed) }
  df <- genGradient(x, func, h=h, diff.type=diff.type)
  return(df)
}

##' @rdname loglikeSTGrad
##' @export
loglikeSTHessian <- function(x, STmodel, type="p", x.fixed=NULL, h=1e-3){
  func <- function(x0){ loglikeST(x0, STmodel, type, x.fixed) }
  H <- genHessian(x, func, h=h)
  return(H)
}

##' @rdname loglikeSTGrad
##' @export
loglikeSTnaiveGrad <- function(x, STmodel, type="p", x.fixed=NULL,
                              h=1e-3, diff.type=0){
  func <- function(x0){ loglikeSTnaive(x0, STmodel, type, x.fixed) }
  df <- genGradient(x, func, h=h, diff.type=diff.type)
  return(df)
}

##' @rdname loglikeSTGrad
##' @export
loglikeSTnaiveHessian <- function(x, STmodel, type="p", x.fixed=NULL, h=1e-3){
  func <- function(x0){ loglikeSTnaive(x0, STmodel, type, x.fixed) }
  H <- genHessian(x, func, h=h)
  return(H)
}
