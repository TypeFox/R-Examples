#'Virtua VECM model
#'
#'Pedagogical tool to create a symbolic VECM model, i.e. just for
#'representation purpose.
#'
#'
#'@param alpha Matrix of alpha speed adjustment coefficients.
#'@param beta Matrix of alpha, cointegrating coefficients.
#'@param lags Matrix containg the lags coefficients.
#'@param inc Matrix containg the include (see following arg) coefficients.
#'@param include Character indicating the type of deterministic term included,
#'if any.
#'@return An object of class \sQuote{VECM}, without however any data.
#'@keywords ts VECM cointegration
#'@export
#'@examples
#'
#'
#'  a<-matrix(c(-0.4, 0.1), ncol=1)
#'  b<-matrix(c(1, -2), ncol=2)
#'
#'  # VECM_symb(alpha=a, beta=t(b))
#'  d<- VECM_symbolic(alpha=a, beta=t(b))
#'  VARrep(d)
#'  d<- VECM_symbolic(alpha=a, beta=t(b), lags=matrix(0, ncol=2, nrow=2))
#'  VARrep(d)
#'  LagMat <- matrix(c(0.1, 0.3, 0.1, 0.2), ncol=2, nrow=2)
#'  incMat <- matrix(c(0.5, 0.1), ncol=1)
#'  d3<- VECM_symbolic(alpha=a, beta=t(b), lags=LagMat, inc=incMat, include="const")
#'  VARrep(d3)
#'
#'
VECM_symbolic <- function(alpha, beta, lags, inc, include = c("none", "const", "trend", "both")){

  include <- match.arg(include)

## few checks
  if(!all(sapply(list(alpha,beta), is.matrix))) stop("args alpha and betas should be matrices")
  if(!all(dim(alpha)==dim(beta))) stop("'alpha' and 'beta' should be same dimension")

## guess values of K and r
  K <- nrow(alpha)
  r <- ncol(alpha)


  if(!missing(lags)){
    if(!is.matrix(lags)) stop("arg 'lags' should be matrix")
    lags_dim <- 1
    lags_mat <- lags
  } else {
      lags_dim <- 0
      lags_mat <- NULL
  }


  if(!missing(inc)){
    if(!is.matrix(inc)) stop("arg 'inc' should be matrix")
    if(nrow(inc)!=K) stop("arg 'inc' should have ",  K, " rows")
    ninc <- switch(include, none=0, const=1, trend=1, both=2)
    inc_nam <- switch(include, none=NULL, const="Intercept", trend="Trend", both=c("Intercept", "Trend"))
    if(ncol(inc)!=ninc) stop("arg 'inc' should have",  K, "cols")
    inc_mat <- inc
  } else {
    inc_mat <- inc_nam <- NULL
  }
  
  Pi <- alpha %*%t(beta)

## build coefs matrices:

  co <- cbind(inc_mat,alpha, lags_mat)
  ECTnam <- if(r==1) "ECT" else paste("ECT", 1:r, sep="")
  lagsnam <- if(lags_dim!=0) paste(paste("X", 1:K, sep=""), rep(1:max(1, lags_dim), each=K), sep=" -") else NULL
  colnames(co) <- c(inc_nam, ECTnam, lagsnam)

## Return results
  res <- list()
  res$lag <- lags_dim
  res$k <- K
  res$model.specific$r <- r
  res$model.specific$beta <- beta
  res$coefficients <- co
  res$include <- include
  res$exogen <- FALSE
  res$model.specific$LRinclude <- "none"

  class(res) <- c("VECM_symb", "VECM")
  return(res)
}


#### TEST
if(FALSE){
  library(tsDyn)

  a<-matrix(c(-0.4, 0.1), ncol=1)
  b<-matrix(c(1, -2), ncol=2)

  # VECM_symbolic(alpha=a, beta=t(b))
  d<- VECM_symbolic(alpha=a, beta=t(b))
  VARrep(d)
  d<- VECM_symbolic(alpha=a, beta=t(b), lags=matrix(0, ncol=2, nrow=2))
  VARrep(d)
  d3<- VECM_symbolic(alpha=a, beta=t(b), lags=matrix(c(0.1, 0.3, 0.1, 0.2), ncol=2, nrow=2), inc=matrix(c(0.5, 0.1), ncol=1), include="const")
  VARrep(d3)



}
