nlpMarginalSkewnorm <- function(sel, y, x, priorCoef=momprior(tau=0.348), priorVar=igprior(alpha=0.01,lambda=0.01), priorSkew=momprior(tau=0.348), method='Laplace', B=10^5, logscale=TRUE, XtX, ytX) {
#Marginal density of the data under non-local priors and skew-normal residuals
# - sel: vector with indexes of variables included in the model
# - y: response variable
# - x: design matrix
  if ((!is.logical(sel)) || (!is.vector(sel))) stop("sel must be a logical vector")
  if (priorCoef@priorDistr != priorSkew@priorDistr) stop("Prior in priorCoef and priorSkew must have same functional form (e.g. mom)") 
  if (priorCoef@priorDistr=='pMOM') {
    r <- as.integer(priorCoef@priorPars['r']); prior <- as.integer(0)
  } else if (priorCoef@priorDistr=='piMOM') {
    r <- as.integer(1); prior <- as.integer(1)
  } else if (priorCoef@priorDistr=='peMOM') {
    r <- as.integer(1); prior <- as.integer(2)
  } else if (priorCoef@priorDistr=='zellner') {
    r <- as.integer(1); prior <- as.integer(3)
  } else {
    stop('Prior specified in priorDistr not recognized')
  }
  tau <- as.double(priorCoef@priorPars['tau'])
  taualpha <- as.double(priorSkew@priorPars['tau'])
  alpha <- as.double(priorVar@priorPars['alpha']); lambda <- as.double(priorVar@priorPars['lambda'])
  #  
  if (is.matrix(y)) y <- as.vector(y)
  if (is.vector(x)) { x <- matrix(x,ncol=1) } else { x <- as.matrix(x) }
  if (missing(XtX)) { XtX <- t(x) %*% x } else { XtX <- as.matrix(XtX) }
  if (missing(ytX)) { ytX <- as.vector(matrix(y,nrow=1) %*% x) } else { ytX <- as.vector(ytX) }
  if (is.logical(sel)) sel <- which(sel)
  if ((length(sel)>0) && ((min(sel)<1) | max(sel)>ncol(x))) stop('Invalid specification of parameter sel')
  sel <- as.integer(sel-1); nsel <- as.integer(length(sel)); 
  p <- as.integer(ncol(x)); n <- as.integer(nrow(x))
  y <- as.double(y); sumy2 <- sum(y^2)
  tau <- as.double(tau); r <- as.integer(r)
  if (method=='auto') method=-1 else if (method=='Laplace') method=0 else if (method=='MC') method=1 else if (method=='plugin') method=2 else stop("Invalid 'method'")
  method <- as.integer(method)
  B <- as.integer(B); logscale <- as.integer(logscale)
  alpha <- as.double(alpha); lambda <- as.double(lambda)
  ans <- .Call("nlpMarginalSkewNormI",sel,nsel,n,p,y,sumy2,x,XtX,ytX,tau,taualpha,r,method,B,logscale,alpha,lambda,prior)
  return(ans);
}
