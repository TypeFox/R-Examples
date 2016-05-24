if(getRversion()>="2.15.1") {
  utils::globalVariables(c("na.omit","coef"))
}

#' @title Function that fits the scglr model
#' @description Calculates the components to predict all the dependent variables.
#' @export scglr
#' @importFrom stats model.matrix model.extract coef cor
#' @param formula an object of class \code{Formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data  a data frame to be modeled.
#' @param family a vector of character of the same length as the number of dependent variables: 
#' "bernoulli", "binomial", "poisson" or "gaussian" is allowed.
#' @param K number of components, default is one.
#' @param size describes the number of trials for the binomial dependent variables. 
#' A (number of statistical units * number of binomial dependent variables) matrix is expected.
#' @param offset used for the poisson dependent variables.
#' A vector or a matrix of size: number of observations * number of Poisson dependent variables is expected.
#' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
#' @param na.action a function which indicates what should happen when the data contain NAs. The default is set to \code{na.omit}.
#' @param crit a list of two elements : maxit and tol, describing respectively the maximum number of iterations and 
#' the tolerance convergence criterion for the Fisher scoring algorithm. Default is set to 50 and 10e-6 respectively. 
#' @param method Regularization criterion type. Object of class "method.SCGLR" 
#' built by \code{\link{methodLPLS}} for PLS-type approach or \code{\link{methodSR}} for Structural Relevance.
#' @return an object of the SCGLR class.
#' @return The function \code{\link{summary}} (i.e., \code{\link{summary.SCGLR}}) can be used to obtain or print a summary of the results.
#' @return The generic accessor functions \code{\link{coef}} can be used to extract various useful features of the value returned by \code{scglr}.
#' @return An object of class "\code{SCGLR}" is a list containing following components:
#' @return \item{u}{matrix of size (number of regressors * number of components), contains the component-loadings, 
#' i.e. the coefficients of the regressors in the linear combination giving each component.}
#' @return \item{comp}{matrix of size (number of statistical units * number of components) having the components as column vectors.}
#' @return \item{compr}{matrix of size (number of statistical units * number of components) having the standardized components as column vectors.}
#' @return \item{gamma}{list of length number of dependant variables. Each element is a matrix of coefficients, standard errors, z-values and p-values.}
#' @return \item{beta}{matrix of size (number of regressors + 1 (intercept) * number of dependent variables), contains the coefficients 
#' of the regression on the original regressors X.}
#' @return \item{lin.pred}{data.frame of size (number of statistical units * number of dependent variables), the fitted linear predictor.}
#' @return \item{xFactors}{data.frame containing the nominal regressors.}
#' @return \item{xNumeric}{data.frame containing the quantitative regressors.}
#' @return \item{inertia}{matrix of size (number of components * 2), contains the percentage and cumulative percentage 
#' of the overall regressors' variance, captured by each component.}
#' @return \item{deviance}{vector of length (number of dependent variables), gives the deviance of each \eqn{y_k}'s GLM on the components.}
#' @references Bry X., Trottier C., Verron T. and Mortier F. (2013) Supervised Component Generalized Linear Regression using a PLS-extension of the Fisher scoring algorithm. \emph{Journal of Multivariate Analysis}, 119, 47-60.
#' @examples \dontrun{
#' library(SCGLR)
#' 
#' # load sample data
#' data(genus)
#' 
#' # get variable names from dataset
#' n <- names(genus)
#' ny <- n[grep("^gen",n)]    # Y <- names that begins with "gen"
#' nx <- n[-grep("^gen",n)]   # X <- remaining names
#' 
#' # remove "geology" and "surface" from nx
#' # as surface is offset and we want to use geology as additional covariate
#' nx <-nx[!nx%in%c("geology","surface")]
#' 
#' # build multivariate formula
#' # we also add "lat*lon" as computed covariate
#' form <- multivariateFormula(ny,c(nx,"I(lat*lon)"),c("geology"))
#' 
#' # define family
#' fam <- rep("poisson",length(ny))
#' 
#' genus.scglr <- scglr(formula=form,data = genus,family=fam, K=4,
#'  offset=genus$surface)
#' 
#' summary(genus.scglr)
#' }
scglr <-  function(formula,data,family,K=1,size=NULL,offset=NULL,subset=NULL,na.action=na.omit,crit=list(),method=methodSR())
{
  #todo family "toto" -> rep("toto",ny)
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula","data","size","offset","subset","na.action"),names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  form <- as.Formula(formula)
  mf$formula <- form
  if(!is.null(size))  size <- as.matrix(size)
  
  if(!is.null(offset)) {
    if(is.vector(offset)) {
      offset <- matrix(offset,nrow(data), sum(family=="poisson"))
    } else {
      offset <- as.matrix(offset)
    }
  } 
  
  mf$size <- size
  mf$offset <- offset
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  crit <- do.call("critConvergence", crit)
  y <- as.matrix(model.part(form,data=mf,lhs=1))
  x <- model.part(form, data=mf, rhs = 1)
  if(length(form)[2]==2){
    AX <- model.part(form, data=mf, rhs = 2)
    namesAx <- names(AX)
    AX <- model.matrix(form,data=mf,rhs=2)[,-1]
    if(is.vector(AX)) {
      AX <- matrix(AX,ncol=1)
      colnames(AX) <- namesAx
    }
  }else{
    AX <- NULL
  }
  vnames <- names(x)
  fTypes <- sapply(x,is.factor)
  if(sum(fTypes)>0){
    xFactors <- x[,fTypes,drop=F]
    colnames(xFactors) <- colnames(x[,fTypes,drop=F])
  }else{
    xFactors <- NULL
  }
  xTypes <- sapply(x,is.numeric)
  if(sum(xTypes)>0){
    xNumeric <- wtScale(x[,xTypes],1/nrow(x))
  }else{
    xNumeric <- NULL
  }
  invsqrtm <- metric(as.data.frame(x))
  x <- model.matrix(form,data=mf)[,-1]
  centerx <- apply(x,2,mean)
  nms <- colnames(x)
  xcr <- scale(x,center=TRUE,scale=FALSE)
  xcr <- xcr%*%invsqrtm
  colnames(xcr) <- nms
  
  ### Controls of  dimension between Y and Size, weights and offsets
  ## number of columns in Y  
  ncy <- ncol(y)
  if(length(family)!=ncy){
    stop("Number of dependent variables and family attributs are different!")
  }
  
  
  if("binomial"%in%family){
    if(is.null(size)){
      stop("Number of trials is unknown for bimomial variables!")
    }else{
      if(ncol(size)!=sum("binomial"==family)){
        stop("Number of trials is different from number of bimomial variables!") 
      }else{
        y[,family%in%"binomial"] <- y[,family%in%"binomial"]/size
      }
    }
  }
  
  if(!is.null(model.extract(mf,"offset"))){
    if(ncol(offset)!=sum("poisson"==family)){
      stop("Number of offset and poisson variables are different!")
    }
  }
  ###compute the K components of scglr 
  size <- model.extract(mf,"size")
  offset <- model.extract(mf,"offset")
  
  kComponent.fit <- kComponents(X=xcr,Y=y,AX=AX,K=K,family=family,size=size,
                                offset=offset,crit=crit,method=method)
  if(is.null(AX)){
    gamma.fit <- multivariateGlm.fit(Y=y,comp=kComponent.fit$comp,family=family,
                                     offset=offset,size=size) 
  }else{    
    gamma.fit <- multivariateGlm.fit(Y=y,comp=cbind(kComponent.fit$comp,AX),family=family,
                                     offset=offset,size=size)       
  }
  gamma.coefs <- sapply(gamma.fit,coef)
  deviance <- sapply(gamma.fit,deviance)
  names(deviance) <- colnames(y)
  colnames(gamma.coefs) <- colnames(y)
  beta.coefs <- f2x(Xcr=xcr,centerx=centerx,invsqrtm=invsqrtm,gamma=gamma.coefs[1:(K+1),],
                    u=kComponent.fit$u,comp=kComponent.fit$comp)
  rownames(beta.coefs) <- c("Intercept",colnames(xcr))
  inertia <- cor(as.matrix(xNumeric), kComponent.fit$compr)
  inertia <- inertia^2
  inertia <- colMeans(inertia)
  names(inertia) <- colnames(kComponent.fit$compr)
  beta <- as.data.frame(rbind(beta.coefs,as.matrix(gamma.coefs[-c(1:(K+1)),])))
  rownames(beta) <- c(rownames(beta.coefs),colnames(AX))
  out <- list(
    call=cl,
    u=as.data.frame(kComponent.fit$u),
    comp=as.data.frame(kComponent.fit$comp),
    compr=as.data.frame(kComponent.fit$compr),
    gamma=lapply(gamma.fit, function(x) summary(x)$coefficients),
    beta=beta, 
    lin.pred= as.data.frame(cbind(1,x)%*%beta.coefs),
    xFactors=as.data.frame(xFactors),
    xNumeric=as.data.frame(xNumeric),
    inertia=inertia,
    deviance=deviance
  )
  class(out) <- "SCGLR"
  return(out)
}





