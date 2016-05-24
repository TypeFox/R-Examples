###################################################################################
#' Constructor of LarsPath class
#'
#' This class stores the results of lars and fusion algorithms.
#'
#' \describe{
#'   \item{nbStep}{Number of steps of the algorithm.}
#'   \item{variable}{List of vector of size "step+1". The i+1-th item contains the index of non-zero coefficients at the i-th step.}
#'   \item{coefficient}{List of vector of size "step+1". The i+1-th item contains the non-zero coefficients at the i-th step.}
#'   \item{l1norm}{Vector of length "step+1", containing the L1-norm of the coefficients at each step.}
#'   \item{lambda}{Vector of length "step+1", containing the lambda at each step.}
#'   \item{dropIndex}{Vector of length "step" containing the index of the dropped variable at the i-th step, 0 means no variable has been dropped at this step.}
#'   \item{addIndex}{Vector of length "step" containing the index of the added variable at the i-th step, 0 means no variable has been added at this step.}
#'   \item{mu}{Intercept.}
#'   \item{meanX}{Mean of columns of X.}
#'	 \item{ignored}{A vector containing index of ignored variables during the algorithm.}
#'   \item{p}{Total number of covariates.}
#'	 \item{fusion}{If TRUE,  results from HDfusion function.}
#'   \item{error}{Error message from lars.}
#' }
#'
#' @aliases LarsPath
#' @name LarsPath-class
#' @rdname LarsPath-class
#' @exportClass LarsPath
#'
#' @seealso \code{\link{HDlars}}
#'
setClass(
  Class="LarsPath",
  representation=representation(
    variable="list",
    coefficient="list",
    l1norm="numeric",
    lambda="numeric",
    dropIndex="list",
    addIndex="list",
    nbStep="numeric",
    mu="numeric",
    meanX="numeric",
    ignored="numeric",
    fusion="logical",
    p="numeric",
    error="character"
  ),
  prototype=prototype(
    variable=list(),
    coefficient=list(),
    l1norm=numeric(0),
    lambda=numeric(0),
    dropIndex=list(),
    addIndex=list(),
    nbStep=numeric(0),
    mu=numeric(0),
    meanX=numeric(0),
    ignored=numeric(0),
    fusion=FALSE,
    p=numeric(0),
    error=character()
  )
)


#' 
#' create a matrix with all estimated coefficients from the output of \code{\link{HDlars}} or \code{\link{EMlasso}} functions. 
#'
#' @title List to sparse matrix conversion
#'
#' @param x a \code{\link{LarsPath}} or \code{EMlasso} object
#' @param row if covariates, covariates are in row
#'
#' @return A sparse matrix containing the values of estimated coefficients for all penalty parameter and all covariates
#' 
#' @export
#' 
#' @seealso \code{\link{HDlars}} \code{\link{EMlasso}}
#' 
listToMatrix = function(x, row = c("covariates","lambda"))
{
  if(!(class(x)%in%c("LarsPath","EMlasso")))
    stop("x must be a LarsPath or EM object")
  row = match.arg(row)
  
  if(class(x)=="EMlasso")
  {
    p = x$p
    coefficient = x$coefficient
    variable = x$variable
    lambda = x$lambda
    l1norm = x$lambda
  }
  else
  {
    p = x@p
    coefficient = x@coefficient
    variable = x@variable
    lambda = x@lambda
    l1norm = x@l1norm
  }
  
  if(row == "covariates")
  {
    bet = Matrix(0, nrow = p, ncol = length(l1norm))
    for(i in 1:length(l1norm))
    {
      bet[variable[[i]], i] = bet[variable[[i]], i] + coefficient[[i]]
    }
    
    return(bet)
  }
  else
  {
    bet = Matrix(0, nrow = length(l1norm), ncol = p)
    for(i in 1:length(l1norm))
    {
      bet[i, variable[[i]]] = bet[i, variable[[i]]] + coefficient[[i]]
      bet[i, variable[[i]]] = bet[i, variable[[i]]] + coefficient[[i]]
    }
    
    return(bet)
  }
  
}


###################################################################################
#' 
#'  plot the path of the lars algorithm.
#'  
#' 
#' 
#' @title plot methods for LarsPath object
#' @param x LarsPath object
#' @param sep.line If TRUE, print vertical dashed line when a variable is added or dropped in the path
#' @param abscissa either "l1norm" or "lambda". If "lambda", regularization parameter is used as abscissa, else l1 norm of the solution is used.
#' @param log.scale If TRUE, use logarithm scale on abscissa
#' @param ... Other plot arguments
#' @docType methods
#' @rdname plot-methods
#' @name plot-methods 
#' @aliases plot,LarsPath-method plot-methods
#'
#' @seealso \code{\link{HDlars}} \code{\link{LarsPath}}
#'
#' @export
setMethod(
  f="plot",
  signature="LarsPath",
  definition=function(x, sep.line = FALSE, abscissa = c("l1norm", "lambda"), log.scale = FALSE, ...)
  {
    #check arguments
    abscissa = match.arg(abscissa)
    if(!is.logical(sep.line))
      stop("sep.line must be a boolean.")
    if(!is.logical(log.scale))
      stop("log.scale must be a boolean.")    
    
    beta = listToMatrix(x, "lambda")
    absciss = x@l1norm
    if(abscissa != "l1norm")
      absciss = c(x@lambda,0)
    logs = ifelse(log.scale, "x", "")
    xlabel = ifelse(abscissa == "l1norm", "l1 norm", "lambda")
    xlabel = paste0(ifelse(log.scale,"Log ", ""), xlabel)
    matplot(absciss, beta, type="l", main = "Lars path", log = logs, xlab = xlabel, ylab = "Coefficients")
  }
)



#' Plot of the coefficients of a step
#'
#' @title Plot of coefficients
#' @param x A LarsPath object.
#' @param step The step at which you want to plot the coefficients.
#' @param ylab Name of the y axis.
#' @param xlab Name of the x axis.
#' @param ... Other plot arguments.
#' @examples 
#' dataset=simul(50,1000,0.4,10,50,matrix(c(0.1,0.8,0.02,0.02),nrow=2))
#' result=HDfusion(dataset$data,dataset$response) 
#' plotCoefficient(result,result@@nbStep) #plot coefficients at the last step
#' @export 
#' 
#' @seealso \code{\link{HDlars}} \code{\link{LarsPath}}
#' 
plotCoefficient=function(x,step,ylab="coefficients",xlab="variables",...)
{
  if(missing(x))
    stop("x is missing.")
  if(missing(step))
    stop("step is missing.")
  if(class(x)!="LarsPath")
    stop("x must be a LarsPath object")
  if(!.is.wholenumber(step))
    stop("step must be a positive integer smaller than x@nbStep")
  if( (step<0) || (step>x@nbStep) )
    stop("step must be a positive integer smaller than x@nbStep")
  
  if(x@fusion)
  {
    index=sort(x@variable[[step+1]],index.return=TRUE)$ix
    if(x@variable[[step+1]][index[1]]!=1)
    { plot(c(1,x@variable[[step+1]][index[1]]-1),rep(0,2),xlim=c(1,x@p),ylim=c(min(0,cumsum(x@coefficient[[step+1]][index])),max(0,cumsum(x@coefficient[[step+1]][index]))),type="l",xlab=xlab,ylab=ylab)
    }else{
      plot(NA,xlim=c(1,x@p),ylim=c(min(cumsum(x@coefficient[[step+1]][index])),max(cumsum(x@coefficient[[step+1]][index]))),xlab=xlab,ylab=ylab)}
    
    a=0
    if(length(x@variable[[step+1]])>1)
      for(i in 1:(length(x@variable[[step+1]])-1))
      {
        a=a+x@coefficient[[step+1]][index[i]]
        lines(c(x@variable[[step+1]][index[i]],x@variable[[step+1]][index[i+1]]-1),rep(a,2))
      }
    
    a=a+x@coefficient[[step+1]][index[length(index)]]
    lines(c(x@variable[[step+1]][index[length(index)]],x@p),rep(a,2))
  }
  else
  {
    plot((1:x@p)[-x@variable[[step+1]]],rep(0,x@p-length(x@coefficient[[step+1]])),ylim=c(min(x@coefficient[[step+1]],0),max(x@coefficient[[step+1]])),xlim=c(1,x@p),xlab=xlab,ylab=ylab,...)
    points(x@variable[[step+1]],x@coefficient[[step+1]],...)    
  }
}



#' Get the vector of coefficients at a given step
#'
#' @title get coefficients at a given step.
#' @param x A LarsPath object.
#' @param step The step at which you want to get the coefficients.
#' @return a vector of size p containing the value of coefficients at the desired step.
#' @examples 
#' dataset=simul(50,1000,0.4,10,50,matrix(c(0.1,0.8,0.02,0.02),nrow=2))
#' result=HDfusion(dataset$data,dataset$response)
#' coefficient=coeff(result,result@@nbStep) #get the coefficients
#' @export 
#' 
#' @seealso \code{\link{HDlars}}  \code{\link{HDfusion}} \code{\link{LarsPath}}
#'
coeff=function(x,step)
{
  if(missing(x))
    stop("x is missing.")
  if(missing(step))
    stop("step is missing.")
  if(class(x)!="LarsPath")
    stop("x must be a LarsPath object")
  if(!.is.wholenumber(step))
    stop("step must be a positive integer smaller than x@step")
  if( (step<0) || (step>x@nbStep) )
    stop("step must be a positive integer smaller than x@step")
  
  beta=rep(0,x@p)
  
  if(x@fusion)
  {
    a=0
    index=sort(x@variable[[step+1]],index.return=TRUE)$ix
    for(i in 1:(length(x@variable[[step+1]])-1))
    {
      a=a+x@coefficient[[step+1]][index[i]]
      beta[x@variable[[step+1]][index[i]]:(x@variable[[step+1]][index[i+1]]-1)]=rep(a, x@variable[[step+1]][index[i+1]]-x@variable[[step+1]][index[i]])
      
    }
    a=a+x@coefficient[[step+1]][index[length(index)]]
    beta[x@variable[[step+1]][index[length(index)]]:x@p]=rep(a,x@p-x@variable[[step+1]][index[length(index)]]+1)
  }
  else
  {
    beta[x@variable[[step+1]]]=x@coefficient[[step+1]]
  }
  
  return(beta)
}


#' Compute coefficients at a given level of penalty
#'
#' @title Compute coefficients
#' @author Quentin Grimonprez
#' @param object a LarsParth object
#' @param index If mode ="norm", index represents the l1-norm of the coefficients with which we want to predict.
#'  If mode="fraction", index represents the ratio (l1-norm of the coefficientswith which we want to predict)/(l1-norm maximal of the LarsPath object).
#'  If mode="lambda", index represents the value of the penalty parameter. If mode="step", index represents the numer of the step at which we want coefficients.
#' @param mode "fraction" or "norm" or "lambda" or "step".
#' @param ... other arguments. Not used
#' @return A vector containing the estimated coefficient for index
#' @method coef LarsPath
#' @examples 
#' dataset=simul(50,10000,0.4,10,50,matrix(c(0.1,0.8,0.02,0.02),nrow=2))
#' result=HDlars(dataset$data[1:40,],dataset$response[1:40])
#' coeff=coef(result,0.3,"fraction")
#' @export
#' 
#' @seealso \code{\link{HDlars}} \code{\link{LarsPath}}
#' 
coef.LarsPath=function(object,index=NULL,mode=c("lambda","step","fraction","norm"),...)
{
  mode <- match.arg(mode)
  
  if(missing(object))
    stop("objectx is missing.")
  if(missing(index) || is.null(index))
    stop("index is missing.")
  if(class(object)!="LarsPath")
    stop("object must be a LarsPath object")
  
  beta=rep(0,object@p)
  
  if(mode=="step")
  {
    if(!.is.wholenumber(index))
      stop("index must be a positive integer smaller than object@step")
    if( (index<0) || (index>object@nbStep) )
      stop("index must be a positive integer smaller than object@step")
    
    if(object@fusion)
    {
      a=0
      ind=sort(object@variable[[index+1]],index.return=TRUE)$ix
      for(i in 1:(length(object@variable[[index+1]])-1))
      {
        a=a+object@coefficient[[index+1]][ind[i]]
        beta[object@variable[[index+1]][ind[i]]:(object@variable[[index+1]][ind[i+1]]-1)]=rep(a, object@variable[[index+1]][ind[i+1]]-object@variable[[index+1]][ind[i]])
        
      }
      a=a+object@coefficient[[index+1]][ind[length(ind)]]
      beta[object@variable[[index+1]][ind[length(ind)]]:object@p]=rep(a,object@p-object@variable[[index+1]][ind[length(ind)]]+1)
    }
    else
    {
      beta[object@variable[[index+1]]]=object@coefficient[[index+1]]
    }
  }
  else
  {
    betatemp=computeCoefficients(object,index,mode)
    beta[betatemp$variable]=betatemp$coefficient
  }
  
  return(beta)
}