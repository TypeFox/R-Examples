#' EM algorithm for lasso penalty
#'
#' @title EM algorithm for lasso penalty
#' @author Quentin Grimonprez, Serge Iovleff
#' @param X the matrix (of size n*p) of the covariates.
#' @param y a vector of length n with the response.
#' @param lambda a sequence of l1 penalty regularization term. If no sequence is provided, the function computes his own sequence.
#' @param maxSteps Maximal number of steps for EM algorithm.
#' @param intercept If TRUE, there is an intercept in the model.
#' @param model "linear" or "logistic"
#' @param burn Number of steps before thresholding some variables to zero.
#' @param threshold Zero tolerance. Coefficients under this value are set to zero.
#' @param eps Epsilon for the convergence of the EM algorithm.
#' @param epsCG Epsilon for the convergence of the conjugate gradient.
#' @return A list containing :
#' \describe{
#'   \item{step}{Vector containing the number of steps of the algorithm for every \code{lambda}.}
#'   \item{variable}{List of vector of the same length as \code{lambda}. The i-th item contains the index of non-zero coefficients for the i-th \code{lambda} value.}
#'   \item{coefficient}{List of vector of the same length as \code{lambda}. The i-th item contains the non-zero coefficients for the i-th \code{lambda} value.}
#'   \item{lambda}{Vector containing the \code{lambda} values.}
#'   \item{mu}{Intercept.}
#' }
#'
#' @examples
#' dataset=simul(50,100,0.4,1,10,matrix(c(0.1,0.9,0.02,0.02),nrow=2))
#' result=EMlasso(dataset$data,dataset$response)
#' # Obtain estimated coefficient in matrix format
#' coefficient = listToMatrix(result)
#' @export
#' 
#' @seealso \code{\link{EMcvlasso}}
#' 
EMlasso <- function(X, y, lambda, maxSteps=1000, intercept=TRUE, model=c("linear", "logistic"), burn=50, threshold=1e-8, eps=1e-5, epsCG=1e-8)
{
  #check arguments
  if(missing(X))
    stop("X is missing.")
  if(missing(y))
    stop("y is missing.")
  .check(X,y,maxSteps,eps,intercept)
  
  ## threshold
  if(!is.double(threshold))
    stop("threshold must be a positive real")
  if(threshold<=0)
    stop("threshold must be a positive real")
  
  ## epsCG
  if(!is.double(epsCG))
    stop("epsCG must be a positive real")
  if(epsCG<=0)
    stop("epsCG must be a positive real")
  
  if(missing(lambda))
    lambda=-1#lambda will be generated in C code
  else
  {
    .check.lambda(lambda)
    lambda=sort(lambda)
    lambda=lambda[lambda>0]
  }
  
  #model
  model = match.arg(model)
  if(model == "logistic")
  {
    # check if y contains 0 and 1
    yb = as.factor(y)
    if(nlevels(yb)!=2)
      stop("In the logistic case, y must contain 0 and 1.")
    
    y = as.numeric(yb)-1
  }
  
  # call em algorithm
  val=list()
  if(model=="linear")
    val=.Call("EMlasso",X,y,lambda,intercept,maxSteps,burn,threshold,eps,epsCG,PACKAGE = "HDPenReg")
  else
    val=.Call("EMlogisticLasso",X,y,lambda,intercept,maxSteps,burn,threshold,eps,epsCG,PACKAGE = "HDPenReg")
  
  val$p = ncol(X)
  
  class(val) = "EMlasso"
  
  return(val)
}


#' EM algorithm for fused-lasso penalty
#'
#' @title EM algorithm for fused-lasso penalty
#' @author Quentin Grimonprez, Serge Iovleff
#' @param X the matrix (of size n*p) of the covariates.
#' @param y a vector of length n with the response.
#' @param lambda1 a positive real. Parameter associated with the lasso penalty.
#' @param lambda2 a positive real. Parameter associated with the fusion penalty.
#' @param maxSteps Maximal number of steps for EM algorithm.
#' @param burn Number of steps before regrouping some variables in segment.
#' @param model "linear" or "logistic"
#' @param intercept If TRUE, there is an intercept in the model.
#' @param eps tolerance for convergence of the EM algorithm.
#' @param eps0 Zero tolerance. Coefficients under this value are set to zero.
#' @param epsCG tolerance for convergence of the conjugate gradient.
#' @return A list containing :
#' \describe{
#'   \item{step}{Vector containing the number of steps of the algorithm for every lambda.}
#'   \item{variable}{List of vector of size "step+1". The i+1-th item contains the index of non-zero coefficients at the i-th step.}
#'   \item{coefficient}{List of vector of size "step+1". The i+1-th item contains the non-zero coefficients at the i-th step.}
#'   \item{lambda}{Vector of length "step+1", containing the lambda at each step.}
#'   \item{mu}{Intercept.}
#' }
#'
#' @examples
#' dataset=simul(50,100,0.4,1,10,matrix(c(0.1,0.9,0.02,0.02),nrow=2))
#' result=EMfusedlasso(dataset$data,dataset$response,1,1)
#' 
#' @export
#' 
#' @seealso \code{\link{EMcvfusedlasso}}
#' 
EMfusedlasso <- function(X, y, lambda1, lambda2, maxSteps=1000, burn=50, intercept=TRUE, model=c("linear","logistic"), eps=1e-5, eps0=1e-8, epsCG=1e-8)
{
  #check arguments
  if(missing(X))
    stop("X is missing.")
  if(missing(y))
    stop("y is missing.")
  .check(X,y,maxSteps,eps,intercept)
  
  ## eps0
  if(!is.double(eps0))
    stop("eps0 must be a positive real")
  if(eps0<=0)
    stop("eps0 must be a positive real")
  
  
  ## epsCG
  if(!is.double(epsCG))
    stop("epSCG must be a positive real")
  if(epsCG<=0)
    stop("epsCG must be a positive real")
  
  ##burn
  if(!.is.wholenumber(burn))
    stop("maxSteps must be a positive integer")
  if( (burn<=0) || (burn>=maxSteps) )
    stop("burn must be a positive integer lower than maxSteps.")
  
  #model
  model = match.arg(model)
  if(model == "logistic")
  {
    # check if y contains 0 and 1
    yb = as.factor(y)
    if(nlevels(yb)!=2)
      stop("In the logistic case, y must contain 0 and 1.")
    
    y = as.numeric(yb)-1
  }
  
  # call EM algorithm
  val=list()
  if(model=="linear")
    val=.Call("EMfusedLasso",X,y,lambda1,lambda2,intercept,maxSteps,burn,eps,eps0,epsCG,PACKAGE = "HDPenReg")
  else
    val=.Call("EMlogisticFusedLasso",X,y,lambda1,lambda2,intercept,maxSteps,burn,eps,eps0,epsCG,PACKAGE = "HDPenReg")
  
  val$p = ncol(X)
  
  class(val) = "EMfusedlasso"
  
  return(val)
}

# check lambda for EM
.check.lambda=function(lambda)
{
  ## lambda: vector of real
  if(!is.numeric(lambda) || !is.vector(lambda))
    stop("lambda must be a vector of positive real.")
  if(length(which(lambda<0))>0)
    stop("lambda must contain only positive number.")
}
