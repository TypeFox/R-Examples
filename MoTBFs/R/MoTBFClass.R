#' Class \code{"motbf"}
#' 
#' Defines an object of class \code{"motbf"} and other basis functions for manipulating 
#' \code{"motbf"} objects.
#' 
#' @name Class-MoTBF
#' @rdname Class-MoTBF
#' @param x Preferably, a list containing, an \code{'mte'} or \code{'mop'} univariate expression
#' and other posibles elements like a \code{"numeric"} vector with the domain of the variable, 
#' the number of iterations needed to solve the optimization problem, among others.
#' Any \R object can be entered, but the utility of this function is not to transform
#' objects of other classes into objects of class \code{"motbf"}.
#' @param class By default is \code{"motbf"}.
#' @param ... Additional arguments, not needed for these methods.
#' @seealso \link{asMOPString} and \link{asMTEString}
#' @examples
#' ## Subclass 'MOP'
#' param <- c(1,2,3,4,5)
#' MOPString <- asMOPString(param)
#' fMOP <- motbf(MOPString)
#' print(fMOP) ## fMOP
#' as.character(fMOP)
#' as.list(fMOP)
#' is(fMOP) ## is.motbf(fMOP)
#' 
#' ## Subclass 'MTE'
#' param <- c(6,7,8,9,10)
#' MTEString <- asMTEString(param)
#' fMTE <- motbf(MTEString)
#' print(fMTE) ## MTE
#' as.character(fMTE)
#' as.list(fMTE)
#' is(fMTE) ## is.motbf(fMTE)


#' @rdname Class-MoTBF
#' @export
motbf <- function(x=0)
{
  if(!is.list(x)) x <- list(Function = x)
  result <- x
  class(result) <- "motbf"
  result
}

#' @rdname Class-MoTBF
#' @export
print.motbf <- function(x, ...) print(x[[1]])

#' @rdname Class-MoTBF
#' @export
as.character.motbf <- function(x, ...) as.character(x[[1]])

#' @rdname Class-MoTBF
#' @export
as.list.motbf <- function(x, ...) as.list(x[[1]])

#' @rdname Class-MoTBF
#'@export
is.motbf <- function(x, class="motbf") is(x, class)

#' Summary of an \code{"motbf"} Object
#' 
#' Summarize an \code{"motbf"} object by describing its main points.
#' 
#' @name summary.motbf
#' @rdname summary.motbf
#' @param object An object of class \code{"motbf"}.
#' @param x An object of class \code{"summary.motbf"}.
#' @param ... further arguments passed to or from other methods.
#' @return The summary of an \code{"motbf"} object. It contains a list of
#' elements with the most important information of the object.
#' @seealso \link{univMoTBF}
#' @examples
#' ## Subclass 'MOP'
#' X <- rnorm(1000)
#' P <- univMoTBF(X, POTENTIAL_TYPE="MOP") ## or POTENTIAL_TYPE="MTE"
#' summary(P)
#' attributes(sP <- summary(P))
#' attributes(sP)
#' sP$Function
#' sP$Subclass
#' sP$Iterations
#' 
#' ## Subclass 'MTE'
#' X <- rnorm(1000)
#' P <- univMoTBF(X, POTENTIAL_TYPE="MTE")
#' summary(P)
#' attributes(sP <- summary(P))
#' attributes(sP)
#' sP$Function
#' sP$Subclass
#' sP$Iterations

#' @rdname summary.motbf
#' @method summary motbf
#' @export
summary.motbf <- function(object, ...)
{
  l <- lapply(1:length(object), function(i) object[[i]])
  names(l) <- attributes(object)$names
  l[[length(l)+1]] <- class(object)
  names(l)[length(l)] <- "Class"
  l[[length(l)+1]] <- coef(object)
  names(l)[length(l)] <- "Coeff"
  class(l) <- "summary.motbf"
  l
}

#' @rdname summary.motbf
#' @export
print.summary.motbf <- function(x, ...)
{
  cat("\n MoTBFs FOR UNIVARIATE DISTRIBUTIONS \n")
  cat("\n Model:"); cat("\n", x$Function, "\n")
  cat("\n Class:", x$Class)
  cat("\n Subclass:", x$Subclass, "\n")
  cat("\n Coefficients:"); cat("\n", x$Coeff, "\n")
  if(!is.null(x$Domain)){ 
    cat("\n Domain:")
    cat("\n (",x$Domain[[1]], ", ",x$Domain[[2]],")", sep="","\n")
  }
  if(!is.null(x$Iterations)&&!is.null(x$Iterations)){
    cat("\n Number of Iterations:", x$Iterations, "\n")
    cat("\n Processing Time:", x$Time, attributes(x$Time)$units, "\n")
  }
  
}



#' Class \code{"jointmotbf"}
#' 
#' Defines an object of class \code{"jointmotbf"} and other basis functions for 
#' manipulating \code{"jointmotbf"} objects.
#' 
#' @name Class-JointMoTBF
#' @rdname Class-JointMoTBF
#' @param x Preferably, a list containing, a joint expression
#' and other posibles elements like a \code{"numeric"} matrix with the domain of the variables, 
#' the dimension of the variables, the number of iterations needed to solve the optimization problem,
#' among others. Any \R object can be entered, but the utility of this function is not to transform
#' objects of other classes into objects of class \code{"jointmotbf"}.
#' @param class By default is \code{"jointmotbf"}.
#' @param ... Additional arguments, not needed for these methods.
#' @seealso \link{jointMoTBF}
#' @examples
#' ## n.parameters is the product of the dimensions
#' dim <- c(3,3)
#' param <- seq(1,prod(dim), by=1)
#' ## Joint Function 
#' f <- list(Parameters=param, Dimensions=dim)
#' jointF <- jointMoTBF(f)
#' 
#' print(jointF) ## jointF
#' as.character(jointF)
#' as.list(jointF)
#' is(jointF)
#' is.jointmotbf(jointF)


#' @rdname Class-JointMoTBF
#' @export
jointmotbf <- function(x = 0)
{
  if(!is.list(x)) x <- list(Function = x)
  result <- x
  class(result) <- "jointmotbf"
  result
}

#' @rdname Class-JointMoTBF
#' @export
print.jointmotbf <- function(x, ...) print(x[[1]])

#' @rdname Class-JointMoTBF
#' @export
as.character.jointmotbf <- function(x, ...) as.character(x[[1]])

#' @rdname Class-JointMoTBF
#' @export
as.list.jointmotbf <- function(x, ...) as.list(x[[1]])

#' @rdname Class-JointMoTBF
#' @export
is.jointmotbf <- function(x, class="jointmotbf") is(x, class)

#' Summary of a \code{"jointmotbf"} Object
#' 
#' Summarize a \code{"jointmotbf"} object by describing its main points.
#' 
#' @name summary.jointmotbf
#' @rdname summary.jointmotbf
#' @param object An object of class \code{"jointmotbf"}.
#' @param x An object of class \code{"summary.jointmotbf"}.
#' @param ... further arguments passed to or from other methods.
#' @return The summary of a \code{"jointmotbf"} object. It contains a list of
#' elements with the most important information of the object.
#' @seealso \link{parametersJointMoTBF} and \link{jointMoTBF}
#' @examples
#' ## 1. EXAMPLE
#' X <- rnorm(100)
#' Y <- rexp(100)
#' data <- data.frame(X, Y)
#' dim <- c(3,4)
#' param <- parametersJointMoTBF(data, dimensions=dim)
#' P <- jointMoTBF(param)
#' summary(P)
#' attributes(sP <- summary(P))
#' attributes(sP)
#' sP$Function
#' sP$Domain
#' sP$Iterations
#' 
#' ##############################################################################
#' ## MORE EXAMPLES #############################################################
#' ##############################################################################
#' ## X <- rnorm(100)
#' ## Y <- rexp(100)
#' ## Z <- rnorm(100, mean=1)
#' ## data <- data.frame(X, Y, Z)
#' ## dim <- c(3,2,4)
#' ## param <- parametersJointMoTBF(data, dimensions=dim)
#' ## P <- jointMoTBF(param)
#' ## summary(P)
#' ## attributes(sP <- summary(P))
#' ## sP$Function
#' ## sP$Domain
#' ## sP$Iterations
#' ##############################################################################
#' ##############################################################################

#' @rdname summary.jointmotbf
#' @method summary jointmotbf
#' @export
summary.jointmotbf <- function(object, ...)
{
  l <- lapply(1:length(object), function(i) object[[i]])
  names(l) <- attributes(object)$names
  l[[length(l)+1]] <- class(object)
  names(l)[length(l)] <- "Class"
  l[[length(l)+1]] <- coef(object)
  names(l)[length(l)] <- "Coeff"
  class(l) <- "summary.jointmotbf"
  l
}

#' @rdname summary.jointmotbf
#' @export
print.summary.jointmotbf <- function(x, ...)
{
  cat("\n MoTBFs FOR MULTIVARIATE DISTRIBUTIONS \n")
  cat("\n Model:"); cat("\n", x$Function, "\n")
  cat("\n Class:", x$Class, "\n")
  cat("\n Coefficients:"); cat("\n", x$Coeff, "\n")
  if(!is.null(x$Domain)){ 
    for(i in 1:ncol(x$Domain)){
      cat("\n Domain ", nVariables(x)[i], ":", sep="")
      cat("\n (",x$Domain[1,i], ", ",x$Domain[2,i],")", sep="")
    }
    cat("\n")
  }
  if(!is.null(x$Iterations)&&!is.null(x$Iterations)){
    cat("\n Number of Iterations:", x$Iterations, "\n")
    cat("\n Processing Time:", x$Time, attributes(x$Time)$units, "\n")
  }
  
}

