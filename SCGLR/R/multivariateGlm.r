#' @title Multivariate generalized linear regression
#' @description \code{multivariateGlm} is used to fit multivariate generalized linear models
#' specified by a symbolic formula together with the distributions of the responses. 
#' This function performs a simple GLM fit for each dependent variable with the associated distribution.
#' @export
#' @importFrom stats model.matrix model.extract 
#' @param formula an object of class \code{Formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data the data frame to be modeled.
#' @param family a vector of character giving the family distribution of each response.
#' @param size a matrix giving the number of trials for each Binomial dependent variable
#' ncol(size) must be equal to the number of Binomial variables.
#' @param offset used for the poisson dependent variables.
#' A vector or a matrix of size: number of observations * number of Poisson dependent variables is expected.
#' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
#' @return the list, each item of which is the glm object associated with each response.
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
#' # remove "geology" and "surface" from nx as surface
#' # is offset and we want to use geology as additional covariate
#' nx <-nx[!nx%in%c("geology","surface")]
#' 
#' # build multivariate formula
#' # we also add "lat*lon" as computed covariate
#' form <- multivariateFormula(ny,c(nx,"I(lat*lon)"),c("geology"))
#' 
#' # split genus dataset
#' sub <- sample(1:nrow(genus),100,replace=FALSE)
#' sub_fit <- (1:nrow(genus))[-sub]
#' 
#' # define family 
#' fam <- rep("poisson",length(ny))
#' 
#' # fit the model
#' genus.scglr <- scglr(formula=form, data=genus, family=fam, K=4, 
#'  offset=genus$surface, subset=sub_fit)
#' 
#' # xnew, the design matrix associated to sub-sample used for prediction
#' # note rhs parameter is introduced to take into account that the 
#' # covariate part of the formula is composed of two differents sets
#' xnew <- model.matrix(form, data=genus[sub,], rhs=1:2)[,-1]
#' 
#' # prediction based on the scglr approch
#' pred.scglr <- multivariatePredictGlm(xnew,family=fam,
#'  beta=genus.scglr$beta, offset=genus$surface[sub])
#' cor.scglr <-diag(cor(pred.scglr,genus[sub,ny])) 
#' plot(cor.scglr, col="red",ylim=c(-1,1))
#' 
#' # prediction based on classical poisson glm
#' genus.glm <- multivariateGlm(formula=form, data=genus, family=fam, 
#'  offset=genus$surface, subset=sub_fit)
#' coefs <- sapply(genus.glm,coef)
#' 
#' # rhs parameter is introduced to take into account that the 
#' # covariate part of the formula is composed of two differents sets
#' pred.glm <- multivariatePredictGlm(xnew,family=fam,beta=coefs,
#'  offset=genus$surface[sub])
#' cor.glm <- diag(cor(pred.glm,genus[sub,ny]))
#' 
#' points(cor.glm, col="blue")
#' }
multivariateGlm<-  function(formula,data,family,size=NULL,offset=NULL,subset=NULL)
{
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula","data","size","offset","subset"),names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  #form <- update(Formula(formula),.~.-1)
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
  y <- as.matrix(model.part(form,data=mf,lhs=1))
  x <- model.part(form, data=mf, rhs = 1)
  if(sum(length(form))==3){
    sx <- model.part(form, data=mf, rhs = 2)
    sx <- model.matrix(form,data=mf,rhs=2)[,-1]
  }else{
    sx <- NULL
  }
  vnames <- names(x)
  fTypes <- sapply(x,is.factor)
  if(sum(fTypes)>0){
    xFactors <- x[,fTypes,drop=FALSE]
    colnames(xFactors) <- colnames(x[,fTypes,drop=FALSE])
  }else{
    xFactors <- NULL
  }
  xTypes <- sapply(x,is.numeric)
  if(sum(xTypes)>0){
    xNumeric <- wtScale(x[,xTypes],1/nrow(x))
  }else{
    xNumeric <- NULL
  }
  #invsqrtm <- metric(as.data.frame(x))
  x <- model.matrix(form,data=mf)[,-1]
  #centerx <- apply(x,2,mean)
  nms <- colnames(x)
  #xcr <- scale(x,center=TRUE,scale=FALSE)
  #xcr <- xcr%*%invsqrtm
  xcr <- x
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
  
  if(is.null(sx)){
    gamma.fit <- multivariateGlm.fit(Y=y,comp=xcr,family=family,
                                     offset=offset,size=size) 
  }else{
    gamma.fit <- multivariateGlm.fit(Y=y,comp=cbind(xcr,sx),family=family,
                                     offset=offset,size=size)     
  }
  
  return(gamma.fit)
}
