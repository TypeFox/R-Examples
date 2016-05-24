if(getRversion()>="2.15.1") {
  utils::globalVariables(c("na.omit","coef"))
}

#' Function that fits and selects the number of component by cross-validation.
#' @export 
#' @importFrom stats model.matrix model.extract coef cor
#' @param formula an object of class "Formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data the data frame to be modeled.
#' @param family a vector of character of length q specifying the distributions of the responses. Bernoulli, binomial, poisson and gaussian are allowed.
#' @param K number of components, default is one.
#' @param nfolds number of folds, default is 5. 
#' Although nfolds can be as large as the sample size (leave-one-out CV), 
#' it is not recommended for large datasets. 
#' @param type loss function to use for cross-validation. 
#' Currently six options are available depending on whether the responses are of the same distribution family.
#' If the responses are all bernoulli distributed, then the prediction performance may be measured
#' through the area under the ROC curve: type = "auc"
#' In any other case one can choose among the following five options ("likelihood","aic","aicc","bic","mspe").
#' @param size specifies the number of trials of the binomial variables included in the model.  A (n*qb) matrix is expected
#' for qb binomial variables.
#' @param offset used for the poisson dependent variables.
#' A vector or a matrix of size: number of observations * number of Poisson dependent variables is expected.
#' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
#' @param na.action a function which indicates what should happen when the data contain NAs. The default is set to the \code{na.omit}.
#' @param crit a list of two elements : maxit and tol, describing respectively the maximum number of iterations and 
#' the tolerance convergence criterion for the Fisher scoring algorithm. Default is set to 50 and 10e-6 respectively. 
#' @param mc.cores max number of cores to use when using parallelization (Not available in windows yet and strongly discouraged if in interactive mode).
#' @param method Regularization criterion type. Object of class "method.SCGLR" 
#' built by \code{\link{methodLPLS}} for PLS-type approach or \code{\link{methodSR}} for Structural Relevance.
#' @return  a matrix containing the criterion values for each response (rows) and each number of components (columns).
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
#' # cross validation
#' genus.cv <- scglrCrossVal(formula=form, data=genus, family=fam, K=12, 
#'  offset=genus$surface)
#' 
#' # find best K
#' mean.crit <- t(apply(genus.cv,1,function(x) x/mean(x)))
#' mean.crit <- apply(mean.crit,2,mean)
#' K.cv <- which.min(mean.crit)-1
#' 
#' #plot(mean.crit, type="l")
#' }
scglrCrossVal <-  function(formula,data,family,K=1,nfolds=5,type="mspe",size=NULL,offset=NULL,subset=NULL,
                           na.action=na.omit,crit=list(), method=methodSR(),mc.cores=1) {
  if( (mc.cores>1) && ((.Platform$OS.type == "windows") || (!requireNamespace("parallel", quietly=TRUE)))){
    warning("Sorry as parallel package is not available, I will use only one core!")
    mc.cores <- 1
  }
  if((mc.cores>1) && interactive()) {
    warning("Using parallel package in interactive mode is strongly discouraged!")
  }
  
  checkLossFunction(type)
  
  if((type=="auc") && (prod(family=="bernoulli")==0))
    stop("auc loss function only when all bernoulli!")
  
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula","data","size","offset","subset","na.action"),names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  #form <- update(Formula(formula),.~.-1)

  form <- as.Formula(formula)
  mf$formula <- form
  if(!is.null(size))  size <- as.matrix(size)
  mf$size <- size
  if(!is.null(offset)) {
    if(is.vector(offset)) {
      offset <- matrix(offset,nrow(data), sum(family=="poisson"))
    } else {
     offset <- as.matrix(offset)
    }
  } 
  
  mf$offset <- offset
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame()) 
  crit <- do.call("critConvergence", crit)
  y <- as.matrix(model.part(form,data=mf,lhs=1))
  x <- model.part(form, data=mf, rhs = 1)
  if(sum(length(form))==3){
    AX <- model.part(form, data=mf, rhs = 2)
    namesAx <- names(AX)
    AX <- model.matrix(form,data=mf,rhs=2)[,-1]
    #browser()
    if(is.vector(AX)) {
      AX <- matrix(AX,ncol=1)
      colnames(AX) <- namesAx
    }
  }else{
    AX <- NULL
  }
  vnames <- names(x)
  invsqrtm <- metric(as.data.frame(x))
  xdesign <- model.matrix(form,data=mf)[,-1]
  centerx <- apply(xdesign,2,mean)
  nms <- colnames(xdesign)
  xcr <- scale(xdesign,center=TRUE,scale=FALSE)
  xcr <- xcr%*%invsqrtm
  colnames(xcr) <- nms
  
  ### Controls of  dimension between Y and Size, weights and offsets
  ## number of columns in Y  
  ncy <- ncol(y)
  if(length(family)!=ncy){
    stop("number of dependent variables and family attributs are different!")
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
  nobs <- nrow(y)
  ny <- ncol(y)
  foldid = sample(rep(seq(nfolds), length = nobs))
  cv <- array(0,c(ny,K,nfolds))
  #logLik <- matrix(0,nfolds,K)
  cvNull <- matrix(0,ny,nfolds)
  
  mainFolds <- function(nf) {
    estid <- which(foldid!=nf)
    valid <- which(foldid==nf)

    cv <- array(0,c(ny,K))
    cvNull <- matrix(0,ny)
    
    kComponent.fit <- kComponents(X=xcr[estid,,drop=FALSE],Y=y[estid,,drop=FALSE],AX=AX[estid,,drop=FALSE],K=K,
                                  family=family,size=size[estid,,drop=FALSE],
                                  offset=offset[estid,,drop=FALSE],crit=crit,method=method)  
    for(kk in seq(K)){ 
      if(is.null(AX)){
        gamma.fit <- multivariateGlm.fit(Y=y[estid,,drop=FALSE],comp=kComponent.fit$comp[,1:kk,drop=FALSE],
                                         family=family,offset=offset[estid,,drop=FALSE],size=size[estid,,drop=FALSE])
      }else{
        gamma.fit <- multivariateGlm.fit(Y=y[estid,,drop=FALSE],comp=cbind(kComponent.fit$comp[,1:kk,drop=FALSE],AX[estid,,drop=FALSE]),
                                         family=family,offset=offset[estid,,drop=FALSE],size=size[estid,,drop=FALSE])  
      }
      gamma.coefs <- sapply(gamma.fit, coef)
      
      beta.coefs <- f2x(Xcr=xcr[estid,,drop=FALSE],centerx=centerx,invsqrtm=invsqrtm,gamma=gamma.coefs[1:(kk+1),,drop=FALSE],u=kComponent.fit$u[,1:kk,drop=FALSE],
                        comp=kComponent.fit$comp[,1:kk,drop=FALSE])
      
      if(is.null(AX)){
        predict <- multivariatePredictGlm(Xnew=x[valid,,drop=FALSE],family,
                                          beta.coefs[,,drop=FALSE],offset[valid,,drop=FALSE])
      }else{
        beta.coefs <- rbind(beta.coefs,gamma.coefs[-c(1:(kk+1)),,drop=FALSE])
        predict <- multivariatePredictGlm(Xnew=cbind(x[valid,,drop=FALSE],AX[valid,,drop=FALSE]),family,
                                          beta.coefs[,,drop=FALSE],offset[valid,,drop=FALSE])
      }
      
      if(type=="auc"){
        cvNull[1:ny] <- auc(roc(y[valid,,drop=FALSE],predict))
      } else if(type%in%c("likelihood","aic","bic","aicc","mspe")){
        cv[1:ny,kk]<- infoCriterion(ynew=y[valid,,drop=FALSE],predict,family,
                                    type=type,size=size[valid,,drop=FALSE],npar=nrow(gamma.coefs))
      }
      
      #      }
      ###Caluclation for NULL model
      if(is.null(AX)){
        gamma.fit <- multivariateGlm.fit(Y=y[estid,,drop=FALSE],comp=NULL,family=family,
                                         offset=offset[estid,,drop=FALSE],
                                         size=size[estid,,drop=FALSE])
        xnew <- rep(1,length(valid))
      }else{
        gamma.fit <- multivariateGlm.fit(Y=y[estid,,drop=FALSE],comp=AX[estid,,drop=FALSE],family=family,
                                         offset=offset[estid,,drop=FALSE],
                                         size=size[estid,,drop=FALSE]) 
        xnew <- AX[valid,,drop=FALSE]
      }
      
      beta.coefs <- sapply(gamma.fit, coef)
      predict <- multivariatePredictGlm(Xnew=xnew,family,
                                        beta.coefs,offset[valid,,drop=FALSE]) 
      
      if(type=="auc"){
        cvNull[1:ny] <- auc(roc(y[valid,,drop=FALSE],predict))
      } else if(type%in%c("likelihood","aic","bic","aicc","mspe")){
        cvNull[1:ny] <- infoCriterion(ynew=y[valid,,drop=FALSE],predict,family,
                                      type=type,size=size[valid,,drop=FALSE],npar=nrow(gamma.coefs))
      }
    }
    return(list(cv=cv,cvNull=cvNull))
  }

  if(mc.cores>1) {
    result <- parallel::mclapply(seq(nfolds),mainFolds,mc.silent=TRUE,mc.cores=mc.cores)
  } else {
    result <- lapply(seq(nfolds),mainFolds)
  }
  cv <- sapply(result,function(x) x$cv,simplify=FALSE)
  cv <- simplify2array(cv)
  cv <- apply(cv,c(1,2),mean)
  
  cvNull <- sapply(result,function(x) x$cvNull,simplify=FALSE)
  cvNull <- simplify2array(cvNull)
  if(is.vector(cvNull)) {
    cvNull <- mean(cvNull)
  } else {
    cvNull <- apply(cvNull,1,mean)
  }
  cv <- cbind(cvNull,cv)
  colnames(cv) <- c("null model",paste("nc",1:K,sep=""))  
  rownames(cv) <- colnames(y)
  return(cv)
}
