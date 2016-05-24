#' liv S4 Object
#'
#' This class is used to store and further analyze the results of the liv function
#' @slot formula  returns an object of type 'formula' as used in the call of the function. Example \code{var1 ~ var2}. 
#' @slot coefficients  model's coefficients
#' @slot seCoefficients  the standard errors of the coefficients
#' @slot groupMeans  the coefficients of the means of the two groups considered to build the latent IV.
#' @slot seMeans  the standard errors of the groups means coefficients.
#' @slot sigma  the coefficients of the variance -covariance matrix.
#' @slot probG1  the coefficient of the probability of group 1.
#' @slot seProbG1  the standard error of the coefcients of probability of group 1.
#' @slot initValues  the initial parameter values.
#' @slot value  the value of the log=likelihood function computed at the optimal parameter values.
#' @slot convCode  the converge code.
#' @slot hessian  the hessian matrix.

#' @name liv-class
#' @rdname liv-class
#' @exportClass liv
#'
#' @examples
#' getSlots("liv")
#'
#' @importFrom methods setClass
#' @export
setClass(

  #Class name
 "liv",


  #Slots / member vars
  slots = c(
    formula = "formula",
    coefficients = "numeric",
    seCoefficients = "numeric",
    groupMeans = "numeric",
    seMeans = "numeric",
    sigma = "matrix",
    probG1 = "numeric",
    seProbG1 = "numeric",
    initValues = "numeric",
    value = "numeric",
    convCode = "integer",
    hessian = "matrix"

  ),

  prototype = list(formula = y ~ x,
                   coefficients = NA_real_,
                   seCoefficients = NA_real_,
                   groupMeans = NA_real_,
                   seMeans = NA_real_,
                   sigma = matrix(NA),
                   probG1 = NA_real_,
                   initValues = NA_real_,
                   value = NA_real_,
                   convCode = NA_integer_,
                   hessian = matrix(NA)
  )
)


#' S3 Method for LIV object for generic "coef"
#'@param object an object of class "liv", usually, a result of a call to liv().
#'@param ... further arguments passed to or from other methods.
#'@export
coef.liv <- function(object, ...)
  {print(object@coefficients)}


#' S3 Method for liv object for generic "summary"
#'@param object an object of class "liv", usually, a result of a call to liv().
#'@param ... further arguments passed to or from other methods.
#'@export
summary.liv <- function(object, ...)
{
    z <- object
    est <- z@coefficients  # estimates value
    se <- z@seCoefficients  # standard errors
    names.coef <- all.vars(z@formula[[3]])

    coef.table <- cbind(est,se)
    colnames(coef.table) <- c("Estimate","Std. Error")
    rownames(coef.table) <- c("Intercept",names.coef)

    cat("\nCoefficients:\n")
    stats::printCoefmat(coef.table) # print the coefficient and std errors

    cat("\nInitial Parameter Values:\n", z@initValues)  # print initial param values
    cat("\n")
    cat("\nThe Value of the log likelihood function:\n", z@value)  # print logLik values
    cat("\n")
    cat("\nConvergence Code:\n", z@convCode)  # print comvergence code
}


#' S3 Method for liv object for generic "print"
#'@param x an object of class "liv", usually, a result of a call to liv().
#'@param ... further arguments passed to or from other methods.
#'@export
print.liv <- function(x, ...)
{
 utils::str(x)
  
}


