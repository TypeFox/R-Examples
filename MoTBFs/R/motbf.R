#' Fitting MoTBFs
#' 
#' Function for fitting univariate mixture of truncated basis functions.
#' Least square optimization is used to minimize the quadratic 
#' error between the empirical cumulative distribution and the estimated one. 
#' 
#' @param data A \code{"numeric"} data vector.
#' @param POTENTIAL_TYPE A \code{"character"} string specifying the potential
#' type, must be one of \code{"MOP"} or \code{"MTE"}.
#' @param evalRange A \code{"numeric"} vector with the range where defining
#' the function. By default: it is \code{NULL} and the function is defined
#' over the data range.
#' @param nparam Number of parameters of the function. By default: it is \code{NULL}
#' and the function fits the best MoTBF taking into account the Bayesian information 
#' criterion (BIC) to penalize the functions. It evaluates the two next functions,
#' if the BIC value doesn't improve then the function with the last best BIC is returned.
#' @param maxParam A \code{"numeric"} value which indicate the maximum number of coefficients in the function. By default it is \code{NULL}; 
#' if not, the output is the function which gets the best BIC with at most this number of parameters.
#' @return An \code{"motbf"} function learned from data.
#' @export
#' @examples
#' ## 1. EXAMPLE
#' ## Data
#' X <- rnorm(5000)
#' 
#' ## Learning
#' f1 <- univMoTBF(X, POTENTIAL_TYPE = "MTE"); f1
#' f2 <- univMoTBF(X, POTENTIAL_TYPE = "MOP"); f2
#' 
#' ## Plots
#' hist(X, prob = TRUE, main = "")
#' plot(f1, xlim = range(X), col = 1, add = TRUE)
#' plot(f2, xlim = range(X), col = 2, add = TRUE)
#' 
#' ## Data test
#' Xtest <- rnorm(1000)
#' ## Filtered data test
#' Xtest <- Xtest[Xtest>=min(X) && Xtest<=max(X)]
#' 
#' ## Log-likelihood
#' sum(log(as.function(f1)(Xtest)))
#' sum(log(as.function(f2)(Xtest)))
#' 
#' ## 2. EXAMPLE
#' ## Data
#' X <- rchisq(5000, df = 5)
#' 
#' ## Learning
#' f1 <- univMoTBF(X, POTENTIAL_TYPE = "MTE", nparam = 11); f1
#' f2 <- univMoTBF(X, POTENTIAL_TYPE = "MOP", maxParam = 10); f2
#' 
#' ## Plots
#' hist(X, prob = TRUE, main = "")
#' plot(f1, xlim = range(X), col = 3, add = TRUE)
#' plot(f2, xlim = range(X), col = 4, add = TRUE)
#' 
#' ## Data test
#' Xtest <- rchisq(1000, df = 5)
#' ## Filtered data test
#' Xtest <- Xtest[Xtest>=min(X) && Xtest<=max(X)]
#' 
#' ## Log-likelihood
#' sum(log(as.function(f1)(Xtest)))
#' sum(log(as.function(f2)(Xtest)))
#' 
univMoTBF <- function(data, POTENTIAL_TYPE, evalRange=NULL, nparam=NULL,  maxParam=NULL)
{
  if(is.null(evalRange)) evalRange <- range(data) else evalRange <- evalRange
  if(is.null(nparam)){
    if(POTENTIAL_TYPE=="MOP"){
      P=bestMOP(data, evalRange, maxParam=maxParam)$bestPx
    } else if(POTENTIAL_TYPE=="MTE"){
      P=bestMTE(data, evalRange, maxParam=maxParam)$bestPx
    } else{
      return(cat("Unknown method, please use MOP or MTE \n"))
    } 
  } else {
    if(POTENTIAL_TYPE=="MOP"){
      P <- mop.learning(data, nparam, evalRange)
    } else if(POTENTIAL_TYPE=="MTE"){
      P <- mte.learning(data, nparam, evalRange)
    } else{
      return(cat("Unknown method, please use MOP or MTE"))
    } 
  }
  return(P)
}


#'Computing the BIC Score of an MoTBF Function
#'
#'Computes the Bayesian information criterion value (BIC) of a 
#'mixture of truncated basis function. The BIC score is the log likelihood 
#'penalizes by the number of parameters of the function and the number of
#'records of the evaluated data.
#'
#'@param Px A function of class \code{"motbf"}.
#'@param X A \code{"numeric"} vector with the data to evaluate.
#'@return A \code{"numeric"} value containing the BIC score.
#'@seealso \link{univMoTBF}
#'@export
#'@examples
#'
#'## Data
#'X <- rexp(10000)
#'
#'## Data test
#'Xtest <- rexp(1000)
#'Xtest <- Xtest[Xtest>=min(X) && Xtest<=max(X)]
#'
#'## Learning
#'f1 <- univMoTBF(X, POTENTIAL_TYPE = "MOP", nparam = 10); f1
#'f2 <- univMoTBF(X, POTENTIAL_TYPE = "MTE", maxParam = 11); f2
#'
#'## BIC values
#'BICMoTBF(Px = f1, X = Xtest)
#'BICMoTBF(Px = f2, X = Xtest)
#'
#'
BICMoTBF <- function(Px, X){
  pPx  <-  as.function(Px)(X)
  size  <-  length(coef(Px))
  BiC  <- sum(log(pPx))-(1/2*(size+1)*log(length(X)))
  return(BiC)
}


#' Extract MoTBF Coefficients
#' 
#' Extracts the parameters of the learned mixtures of truncated basis
#' functions. 
#' 
#' @param object An MoTBF function.
#' @param \dots other arguments
#' @return A numeric vector with the parameters of the function.
#' @seealso \link{univMoTBF}, \link{coeffMOP} and \link{coeffMTE}
#' @export
#' @examples
#'
#'## Data
#'X <- rchisq(2000, df = 5)
#'
#'## Learning
#'f1 <- univMoTBF(X, POTENTIAL_TYPE = "MOP"); f1
#'## Coefficients
#'coef(f1)
#'
#'## Learning
#'f2 <- univMoTBF(X, POTENTIAL_TYPE = "MTE", maxParam = 10); f2
#'## Coefficients
#'coef(f2)
#'
#'## Learning
#'f3 <- univMoTBF(X, POTENTIAL_TYPE = "MOP", nparam=10); f3
#'## Coefficients
#'coef(f3)
#'
#'## Plots
#'plot(NULL, xlim = range(X), ylim = c(0,0.2), xlab="X", ylab="density")
#'plot(f1, xlim = range(X), col = 1, add = TRUE)
#'plot(f2, xlim = range(X), col = 2, add = TRUE)
#'plot(f3, xlim = range(X), col = 3, add = TRUE)
#'hist(X, prob = TRUE, add= TRUE)
#'
coef.motbf <- function(object, ...)
{
  if(is.mop(object)) parameters <- coeffMOP(object)
  if(is.mte(object)) parameters <- coeffMTE(object)
  return(parameters)
}


#' Subclass \code{"motbf"} Functions
#'
#' Collection of functions for detecting the subclass of an \code{"motbf"}
#' object. It can be \code{"mop"} or \code{"mte"}.
#' 
#' @name Subclass-MoTBF
#' @rdname Subclass-MoTBF
#' @param fx A function of the class \code{"motbf"}.
#' @return \code{is.mte} and \code{is.mop} return a logical value, \code{TRUE} if it is an \code{"motbf"} object of the subclass
#' \code{"mte"} or \code{"mop"}, respectly; or \code{FALSE} otherwise. 
#' \code{subclass} returns a \code{"character"} string, \code{"mte"} or \code{"mop"}.
#' @seealso \link{univMoTBF}
#' @examples
#' 
#' ## MOP Function
#' X <- rnorm(1000)
#' P <- univMoTBF(X, POTENTIAL_TYPE="MOP")
#' is.mop(P)
#' subclass(P)
#' 
#' ## MTE Function
#' X <- rchisq(1000, df=4)
#' P <- univMoTBF(X, POTENTIAL_TYPE="MTE")
#' is.mte(P)
#' subclass(P)

#' @rdname Subclass-MoTBF
#' @export
is.mte <- function(fx)
{
  if(!is.null(fx$Subclass)) return(fx$Subclass=="mte")
  f <- fx[[1]]; l <- length(strsplit(f, split="exp", fixed=T)[[1]])-1
  return(l!=0&&is.motbf(fx))
}

#' @rdname Subclass-MoTBF
#' @export
is.mop <- function(fx)
{
  if(!is.null(fx$Subclass)) return(fx$Subclass=="mop")
  f <- fx[[1]]; l <- length(strsplit(f, split="exp", fixed=T)[[1]])-1
  return(l==0&&is.motbf(fx))
}

#' @rdname Subclass-MoTBF
#' @export
subclass <- function(fx)
{
  if(is.mop(fx)) return("mop")
  if(is.mte(fx)) return("mte")
}


#' Derivative MoTBF
#' 
#' Compute derivatives of \code{"motbf"} objects.
#' 
#' @param fx An object of \code{"motbf"} class.
#' @return The derivative which is also an \code{"motbf"} function.
#' @seealso \link{univMoTBF}, \link{derivMOP} and \link{derivMTE}
#' @export
#' @examples
#' 
#' ## 1. EXAMPLE
#' X <- rexp(1000)
#' Px <- univMoTBF(X, POTENTIAL_TYPE="MOP")
#' derivMoTBF(Px)
#' 
#' ## 2. EXAMPLE
#' X <- rnorm(1000)
#' Px <- univMoTBF(X, POTENTIAL_TYPE="MOP")
#' derivMoTBF(Px)
#' 
#' ## 3. EXAMPLE
#' X <- rchisq(1000, df = 3)
#' Px <- univMoTBF(X, POTENTIAL_TYPE="MTE")
#' derivMoTBF(Px)
#' 
#' ## 4. EXAMPLE
#' Px <- "x+2"
#' class(Px)
#' derivMoTBF(Px)
#' ## Message: "It is not an 'motbf' function."
#' 
derivMoTBF <- function(fx)
{
  if(!is.motbf(fx)) return(cat("It is not an 'motbf' function."))
  if(is.mop(fx)) f <- derivMOP(fx)
  if(is.mte(fx)) f <- derivMTE(fx)
  return(f)
}


#' Integral MoTBF
#' 
#' Gets the integral of a one dimensional mixture of truncated basis function 
#' over a bounded or unbounded interval.
#' 
#' @param fx An \code{"motbf"} function.
#' @param min The lower limit of integration. By default it is NULL.
#' @param max The higther limit of the interval. By default it is NULL.
#' @details If the limits of the interval, min and max are NULL, then the output is
#' the expression of the non-defined integral. If only min contains a numeric value,
#' then the expression of the integral is evaluated in this point.
#' @return The non-defined integral of the MoTBF function that is also an \code{"motbf"}
#' function or the defined integral that is a \code{"numeric"} value.
#' @seealso \link{univMoTBF}, \link{integralMOP} and \link{integralMTE}
#' @export
#' @examples
#' 
#' ## 1. EXAMPLE
#' X <- rexp(1000)
#' Px <- univMoTBF(X, POTENTIAL_TYPE="MOP")
#' integralMoTBF(Px)
#' integralMoTBF(Px, 1.2)
#' integralMoTBF(Px, min(X), max(X))
#' 
#' ## 2. EXAMPLE
#' X <- rnorm(1000)
#' Px <- univMoTBF(X, POTENTIAL_TYPE="MOP")
#' iP <- integralMoTBF(Px); iP
#' plot(iP, xlim=range(X))
#' integralMoTBF(Px, 0.2)
#' integralMoTBF(Px, min(X), max(X))
#' 
#' ## 3. EXAMPLE
#' X <- rchisq(1000, df = 3)
#' Px <- univMoTBF(X, POTENTIAL_TYPE="MTE")
#' integralMoTBF(Px)
#' integralMoTBF(Px, 1)
#' integralMoTBF(Px, min(X), max(X))
#' 
#' ## 4. EXAMPLE
#' Px <- "1+x+5"
#' class(Px)
#' integralMoTBF(Px)
#' ## Message: "It is not an 'motbf' function."
#'
integralMoTBF <- function(fx, min=NULL, max=NULL)
{
  if(!is.motbf(fx)) return(cat("It is not an 'motbf' function."))
  if(is.null(min)&&is.null(max)){
    if(is.mop(fx)) return(integralMOP(fx))
    if(is.mte(fx)) return(integralMTE(fx))
  } else if(!is.null(min)&&is.null(max)){
    if(is.mop(fx)) return(as.function(integralMOP(fx))(min))
    if(is.mte(fx)) return(as.function(integralMTE(fx))(min))
  } else{
    if(is.mop(fx)) return(as.function(integralMOP(fx))(max) - as.function(integralMOP(fx))(min))
    if(is.mte(fx)) return(as.function(integralMTE(fx))(max) - as.function(integralMTE(fx))(min))
  } 
}

#' Coerce an \code{"motbf"} Object to a Function
#'
#' Takes an \code{"motbf"} object and contructs an \R function to evaluate it at points.
#' 
#' @param x An object of class \code{"motbf"}.
#' @param ... Further arguments to be passed to or from the method. Not necessary for this method.
#' @details This is an \code{S3} method for the generic function \link{as.function}.
#' @return It returns a function to evaluate an object of class \code{"motbf"}.
#' @export
#' @examples
#' 
#' ## Data
#' X <- rchisq(5000, df = 3)
#' 
#' ## Learning
#' P <- univMoTBF(X, POTENTIAL_TYPE = "MOP"); P
#' 
#' ## Evaluation
#' as.function(P)(min(X))
#' as.function(P)(max(X))
#' as.function(P)(10)
#' density <- as.function(P)(X)
#' 
#' ## Plot
#' hist(X, prob = TRUE, main = "")
#' points(X, density, col=4, pch=16)
#' 
as.function.motbf <- function(x, ...)
{
  if(is.mte(x)){
    string <- x[[1]]
    return(as.function(alist(x = , eval(parse(text = string)))))
  }else{
    v <- nVariables(x)
    string <- x[[1]]
    f <- function(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,
                  s,t,u,v,w,x,y,z) eval(parse(text = string))
    formals(f) <- formals(f)[1:length(v)]
    names(formals(f)) <- v
    return(f)
  } 
}

#' Plots for \code{'motbf'} Objects
#'
#' Draws an \code{'motbf'} function.
#' 
#' @param x An object of class \code{'motbf'}.
#' @param xlim The range to be encompassed by the x axis; by default \code{0:1}.
#' @param ylim The range
#' @param type As for \link{plot}.
#' @param \dots Further arguments to be passed as for \link{plot}.
#' @method plot motbf
#' @return A plot of the specificated function.
#' @export
#' @examples
#'
#'## 1. EXAMPLE 
#'## Data
#'X <- rexp(2000)
#'
#'## Learning
#'f1 <- univMoTBF(X, POTENTIAL_TYPE = "MOP"); f1
#'f2 <- univMoTBF(X, POTENTIAL_TYPE = "MTE", maxParam = 10); f2
#'f3 <- univMoTBF(X, POTENTIAL_TYPE = "MOP", nparam=10); f3

#'## Plots
#'plot(NULL, xlim = range(X), ylim = c(0,0.8), xlab="X", ylab="density")
#'plot(f1, xlim = range(X), col = 1, add = TRUE)
#'plot(f2, xlim = range(X), col = 2, add = TRUE)
#'plot(f3, xlim = range(X), col = 3, add = TRUE)
#'hist(X, prob = TRUE, add= TRUE)
#' 
#'## 2. EXAMPLE 
#'## Data
#'X <- c(rnorm(2000, mean = -3),rnorm(2000, mean = 3))
#'
#'## Learning
#'f1 <- univMoTBF(X, POTENTIAL_TYPE = "MOP"); f1
#'f2 <- univMoTBF(X, POTENTIAL_TYPE = "MTE"); f2

#'## Plots
#'plot(NULL, xlim = range(X), ylim = c(0,0.20), xlab="X", ylab="density")
#'plot(f1, xlim = range(X), col = 2, add = TRUE)
#'plot(f2, xlim = range(X), col = 4, add = TRUE)
#'hist(X, prob = TRUE, add= TRUE)
#' 
plot.motbf <- function(x, xlim=0:1, ylim=NULL, type="l", ...){
  plot(as.function(x), xlim = xlim, ylim = ylim, type="l", ylab="Px", ...)
}


