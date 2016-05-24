sciopath <- function (S, lambdalist=NULL, thr = 1e-4, maxit = 1e4, pen.diag=F, sym=T) {
  p = nrow(S)
  
  if(is.null(lambdalist)){
    lambdalist=seq(max(abs(S))/10,max(abs(S)),length=10)
  }
  lambdalist=sort(lambdalist)
  nlambda <- length(lambdalist)
  
  jerrlist=matrix(0,nlambda, 1)
  nniter <- 1
  ## w <- matrix(0.0,p,p)
  w <- diag(1/diag(S),p)
  wlist <- array(0.0, c(p,p,nlambda))

  idiag <- pen.diag*1
  isym <- sym*1
  
  mode(lambdalist) = "double"
  mode(nlambda) = "integer"
  mode(S) = "double"
  mode(wlist) = "double"
  mode(w) = "double"
  mode(p) = "integer"
  mode(maxit) = "integer"
  mode(nniter) = "integer"
  mode(thr) = "double"
  mode(jerrlist) = "integer"
  mode(idiag) = "integer"
  mode(isym) = "integer"
  
  junk <- .Fortran("sciopath",
                   wlist=wlist,
                   p,
                   S,
                   w,
                   lambdalist,
                   nlambda,
                   thr,
                   maxit,
                   nniter,
                   jerrlist,
                   idiag,
                   isym,
                   PACKAGE="scio")
  wlist <- array(junk$wlist,  c(p,p,nlambda))
  return(list(wlist=wlist, lambdalist=lambdalist))
}

