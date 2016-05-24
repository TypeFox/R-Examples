#' cross validation function for lars algorithm
#'
#' @title cross validation
#' @author Quentin Grimonprez
#' @param X the matrix (of size n*p) of the covariates.
#' @param y a vector of length n with the response.
#' @param nbFolds the number of folds for the cross-validation.
#' @param index Values at which prediction error should be computed. When mode = "fraction", this is the fraction of the saturated |beta|. 
#' The default value is seq(0,1,by=0.01). When mode="lambda", this is values of lambda.
#' @param mode Either "fraction" or "lambda". Type of values containing in partition.
#' @param maxSteps Maximal number of steps for lars algorithm.
#' @param partition partition in nbFolds folds of y. Must be a vector of same size than y containing the index of folds.
#' @param intercept If TRUE, there is an intercept in the model.
#' @param eps Tolerance of the algorithm.
#' @return A list containing 
#' \describe{
#'   \item{cv}{Mean prediction error for each value of index.}
#'   \item{cvError}{Standard error of cv.}
#'   \item{minCv}{Minimal cv criterion.}
#'   \item{minIndex}{Value of index for which the cv criterion is minimal.}
#'   \item{index}{Values at which prediction error should be computed. This is the fraction of the saturated |beta|. The default value is seq(0,1,by=0.01).}
#'   \item{maxSteps}{Maximum number of steps of the lars algorithm.}
#' }
#' @examples 
#' dataset=simul(50,10000,0.4,10,50,matrix(c(0.1,0.8,0.02,0.02),nrow=2))
#' result=HDcvlars(dataset$data,dataset$response,5)
#' @export
#' 
HDcvlars <- function(X,y,nbFolds=10,index=seq(0,1,by=0.01),mode=c("fraction","lambda"),maxSteps=3*min(dim(X)),partition=NULL,intercept=TRUE,eps=.Machine$double.eps^0.5)
{
  #check arguments
  mode <- match.arg(mode)
  if(missing(X))
    stop("X is missing.")
  if(missing(y))
    stop("y is missing.")
  index=unique(index)
  .checkcvlars(X,y,maxSteps,eps,nbFolds,index,intercept,mode)
  
  if(!is.null(partition))
  {
    if(!is.numeric(partition) || !is.vector(partition))
      stop("partition must be a vector of integer.")
    if(length(partition)!=length(y))
      stop("partition and y must have the same size.")
    
    part=table(partition)
    nbFolds=length(part)
    
    if(max(part)-min(part)>1)
      stop("Size of different folds are not good.")
    
    nam=as.numeric(names(part))
    for(i in 1:length(nam))
    {
      if(!(nam[i]==i))
        stop("check the number in the partition vector.")  
    }
    
    #reorder the  partition in decreasing order of size
    ord=order(part,decreasing=TRUE)
    partb=partition
    for(i in 1:nbFolds)
      partition[partb==ord[i]]=i
    
    partition=partition-1
  }
  else
    partition=-1
  
  lambdaMode=FALSE
  if(mode=="lambda")
    lambdaMode=TRUE
  
  # call lars algorithm
  val=.Call( "cvlars",X,y,nrow(X),ncol(X),maxSteps,intercept,eps,nbFolds,partition,index,lambdaMode,PACKAGE = "HDPenReg" )
  
  #create the output object
  cv=list(cv=val$cv,cvError=val$cvError,minCv=min(val$cv),minIndex=index[which.min(val$cv)],index=index,maxSteps=maxSteps,mode=mode)
  
  class(cv)="HDcvlars"
  
  return(cv)
}


# check arguments from cvlars 
.checkcvlars=function(X,y,maxSteps,eps,nbFolds,index,intercept,mode)
{
  ## X: matrix of real
  if(!is.numeric(X) || !is.matrix(X))
    stop("X must be a matrix of real")
  
  ## y: vector of real
  if(!is.numeric(y) || !is.vector(y))
    stop("y must be a vector of real")
  if(length(y)!=nrow(X))
    stop("The number of rows of X doesn't match with the length of y")
  
  ## maxSteps
  if(!.is.wholenumber(maxSteps))
    stop("maxSteps must be a positive integer")
  if(maxSteps<=0)
    stop("maxSteps must be a positive integer")
  
  ## nbFolds
  if(!.is.wholenumber(nbFolds))
    stop("nbFolds must be a positive integer")
  if(nbFolds<=0)
    stop("nbFolds must be a positive integer")
  if(nbFolds>length(y))
    stop("nbFolds must be lower than the number of samples")
  
  ## eps
  if(!is.double(eps))
    stop("eps must be a positive real")  
  if(eps<=0)
    stop("eps must be a positive real")	
  
  ## index
  if(!is.numeric(index) || !is.vector(index))
    stop("index must be a vector")
  if( (mode=="fraction") && (max(index)>1 || min(index)<0))
    stop("index must be a vector of real between 0 and 1")
  if( (mode=="lambda") && (min(index)<0))
    stop("index must be a vector of positive real")
  
  ##intercept
  if(!is.logical(intercept))
    stop("intercept must be a boolean") 
}

#' plot cross validation mean square error
#'
#' @title plot cross validation mean square error
#' @author Quentin Grimonprez
#' @param x Output from HDcvlars function.
#' @param ... graphical parameters
#' @aliases plot.HDcvlars
#' @method plot HDcvlars
#' @examples 
#' dataset=simul(50,10000,0.4,10,50,matrix(c(0.1,0.8,0.02,0.02),nrow=2))
#' result=HDcvlars(dataset$data,dataset$response,5)
#' plot(result)
#' @export
#' 
plot.HDcvlars=function(x,...)
{
  if(missing(x))
    stop("x is missing.")
  if(class(x)!="HDcvlars")
    stop("x must be an output of the HDcvlars function.")
  
  index=x$index
  minIndex=x$minIndex
  lab="Fraction L1 Norm"
  if(x$mode=="lambda")
  {
    lab="log(lambda)"
    index=log(index)
    minIndex=log(minIndex)
  }
  plot(index, x$cv, type = "b", ylim = range(x$cv, x$cv + x$cvError, x$cv - x$cvError),xlab=lab,ylab="Cross-Validated MSE",...)
  lines(index, x$cv+x$cvError,lty=2)
  lines(index, x$cv-x$cvError,lty=2)
  abline(v=minIndex,lty="dotted",col="blue")
  invisible()
}