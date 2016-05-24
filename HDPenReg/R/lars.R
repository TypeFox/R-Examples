#' It performs the lars algorithm for solving lasso problem. 
#' It is a linear regression problem with a l1-penalty on the estimated coefficient. 
#'
#' @title Lars algorithm
#' @author Quentin Grimonprez
#' @param X the matrix (of size n*p) of the covariates.
#' @param y a vector of length n with the response.
#' @param maxSteps Maximal number of steps for lars algorithm.
#' @param intercept If TRUE, add an intercept to the model.
#' @param eps Tolerance of the algorithm.
#' @return An object of type \code{\link{LarsPath}}.
#' @examples 
#' dataset=simul(50,10000,0.4,10,50,matrix(c(0.1,0.8,0.02,0.02),nrow=2))
#' result=HDlars(dataset$data,dataset$response)
#' # Obtain estimated coefficient in matrix format
#' coefficient = listToMatrix(result)
#' 
#' @export
#' 
#' @details
#' The l1 penalty performs variable selection via shrinkage of the estimated coefficient. 
#' It depends on a penalty parameter called lambda controlling the amount of regularization.
#' The objective function of lasso is : \deqn{||y-X\beta||_2 + \lambda||\beta||_1}
#' 
#' @references Efron, Hastie, Johnstone and Tibshirani (2003) "Least Angle Regression" (with discussion) Annals of Statistics
#' 
#' @seealso \code{\link{LarsPath}} \code{\link{HDcvlars}} \code{\link{listToMatrix}}
#' 
HDlars <- function(X,y,maxSteps=3*min(dim(X)),intercept=TRUE,eps=.Machine$double.eps^0.5)
{
  #check arguments
  if(missing(X))
    stop("X is missing.")
  if(missing(y))
    stop("y is missing.")
  .check(X,y,maxSteps,eps,intercept)
  
  # call lars algorithm
  val=.Call( "lars",X,y,nrow(X),ncol(X),maxSteps,intercept,eps,PACKAGE = "HDPenReg" )
  
  #create the output object
  path=new("LarsPath",variable=val$varIdx,coefficient=val$varCoeff,lambda=val$lambda,l1norm=val$l1norm,addIndex=val$evoAddIdx,dropIndex=val$evoDropIdx,
           nbStep=val$step,mu=val$mu,ignored=val$ignored,p=ncol(X),error=val$error,meanX=val$muX)
  return(path)
}

#' It performs the lars algorithm for solving a special case of lasso problem. 
#' It is a linear regression problem with a l1-penalty on the difference of two successive coefficients.
#'
#' @title Fusion algorithm
#' @author Quentin Grimonprez
#' @param X the matrix (of size n*p) of the covariates.
#' @param y a vector of length n with the response.
#' @param maxSteps Maximal number of steps for lars algorithm.
#' @param intercept If TRUE, there is an intercept in the model.
#' @param eps Tolerance of the algorithm.
#' @return An object of type \code{\link{LarsPath}}. \code{\link{LarsPath-class}}.
#' @examples  
#' dataset=simul(50,10000,0.4,10,50,matrix(c(0.1,0.8,0.02,0.02),nrow=2))
#' result=HDfusion(dataset$data,dataset$response)
#' @export
#' 
#' @references Efron, Hastie, Johnstone and Tibshirani (2003) "Least Angle Regression" (with discussion) Annals of Statistics
#' 
#' @seealso LarsPath HDlars
#' 
HDfusion <- function(X,y,maxSteps=3*min(dim(X)),intercept=TRUE,eps=.Machine$double.eps^0.5)
{
  #check arguments
  if(missing(X))
    stop("X is missing.")
  if(missing(y))
    stop("y is missing.")
  .check(X,y,maxSteps,eps,intercept)
  
  # call fusion algorithm
  val=.Call( "fusion",X,y,nrow(X),ncol(X),maxSteps,intercept,eps,PACKAGE = "HDPenReg" )
  
  #create the output object
  path=new("LarsPath",nbStep=val$step,variable=val$varIdx,coefficient=val$varCoeff,lambda=val$lambda,l1norm=val$l1norm,addIndex=val$evoAddIdx,
           dropIndex=val$evoDropIdx,p=ncol(X),fusion=TRUE,error=val$error)
  
  return(path)
}

# check arguments from lars and fusion algorithm
.check=function(X,y,maxSteps,eps,intercept)
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
  
  ## eps
  if(!is.double(eps))
    stop("eps must be a positive real")  
  if(eps<=0)
    stop("eps must be a positive real")	
  
  ##intercept
  if(!is.logical(intercept))
    stop("intercept must be a boolean") 
}

#check if a number is an integer
.is.wholenumber=function(x, tol = .Machine$double.eps^0.5)  
{
  abs(x - round(x)) < tol
}