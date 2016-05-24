# @title Estimation algorithm for K components
# @description calculates the K components by iteratively calling function oneComponent 
# @export
# @param X matrix (n*p) containing the standardized covariates
# @param Y matrix (n*q) containing dependent variables
# @param AX matrix of additional covariates used in the generalized regression but not entering the linear 
# combinations giving components
# @param K integer specifying the number of components
# @param family a vector of the same length as the number of responses containing characters 
# identifying the distribution families of the dependent variables.
# "bernoulli", "binomial", "poisson" or "gaussian" are allowed.
# @param size matrix of size statistical units * number of binomial responses, giving the number of trials 
# for binomial dependent variables. 
# @param offset used for the poisson dependent variables.
# A vector or a matrix of size: number of observations * number of Poisson dependent variables is expected
# @param crit a list of maxit and tol, default is 50 and 10e-6. If responses are bernoulli variables only, tol should generally be increased
# @param method optimization algorithm used for loadings and components. Object of class "method.SCGLR" 
# built by \code{\link{methodEigen}} or \code{\link{methodIng}}
# @return a list 
# @return \item{u}{matrix of size (number of regressors * number of components), contains the component-loadings, 
# i.e. the coefficients of the regressors in the linear combination giving each component}
# @return \item{comp}{matrix of size (number of statistical units * number of components) having the components as column vectors}
# @return \item{compr}{matrix of size (number of statistical units * number of components) having the standardized components as column vectors}
# @return \item{ds}{the final value of the regularization degree}
kComponents <- function(X,Y,AX,K,family,size=NULL,offset=NULL,crit,method)
{
  n<-dim(X)[1]
  p<-dim(X)[2]
  q<-dim(Y)[2]
  if (is.null(q)) q<-1
  ds <- rep(0,K)
  # on stocke la matrice X et on la centre-reduit
  Xorig<-X
  # TODO verbose option
  if(FALSE) print("first component")
  out <- oneComponent(X,Y,AX,family=family,size=size,
                      offset=offset,ds=ds[1],crit=crit,method=method)
  
  if(method$method=="lpls") {
  while (is.logical(out))
  {
    ds[1] <- ds[1]+1
    out <- oneComponent(X,Y,AX,family=family,size=size,
                          offset=offset,ds=ds[1],crit=crit,method=method)
    }	
  }	
  u <- out
  f <- X%*%u
  comp <- f 
  compr <- f/sqrt(c(crossprod(f,f))/n)
  # projection sur f-ortho de X
  Xcour <- X-tcrossprod(f,f)%*%X/c(crossprod(f,f))
  #Composantes suivantes si K>1	
  if(K>1)
  {
    keep_tol <- crit$tol
    for (k in 2:K)
    {
      crit$tol <- keep_tol
      # TODO verbose option
      if(FALSE) message("\nComponent ",k)
      Tab <- cbind(comp,AX)
      out <- oneComponent(Xcour,Y,Tab,family=family,size=size,
                          offset=offset,ds=ds[k],crit=crit,method=method)
      if(method$method=="lpls") {
      while (is.logical(out))
      {
        ds[k] <- ds[k]+1
        if(ds[k]>10){
          ds[k] <- 20+ds[k]
          crit$tol <- min(crit$tol*10,1)
        }
        out <- oneComponent(Xcour,Y,Tab,family=family,size=size,
                              offset=offset,ds=ds[k],crit=crit,method=method)
        }
      }
      if(is.logical(out)) {
        stop("Pas de convergence !!!!!")
      }
      f <- Xcour%*%out
      compr <- cbind(compr,f/sqrt(c(crossprod(f,f))/n))
      u <- cbind(u,out)
      comp <- cbind(comp,f)   
      Xcour <- Xcour-tcrossprod(f,f)%*%Xcour/c(crossprod(f,f))
    }
  }	
  u <- as.matrix(u)
  colnames(u) <- paste("u",1:ncol(u),sep="")
  colnames(comp) <- paste("sc",1:ncol(comp),sep="")
  colnames(compr) <- paste("sc",1:ncol(compr),sep="")

  return(list(u=u,comp=comp,compr=compr,ds=ds))
}


