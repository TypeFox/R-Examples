#' Extract Coefficients and (Co)Variance Matrix from dosresmeta Objects
#'
#' @param object an object of class "\code{dosresmeta}"
#' @param format format of the returned object
#' @param \dots further arguments passed to or from other methods.
#' @return For \code{coef}, a vector (default) or matrix with the estimated (fixed-effects) coefficients.
#' 
#' 
#' For \code{vcov}, the (co)variance matrix of the estimated (fixed-effects) coefficients.
#' 
#' @author Alessio Crippa, \email{alessio.crippa@@ki.se}
#' @seealso \code{\link{dosresmeta}}, \code{\link{coef}}, \code{\link{vcov}}
#' 
#' @examples
#' ## Load data and run the model
#' data("alcohol_cvd")
#' model <- dosresmeta(formula = logrr ~ dose + I(dose^2), type = type, id = id,
#'                    se = se, cases = cases, n = n, data = alcohol_cvd) 
#'
#'## Fixed-effect coefficients
#'coef(model)
#'
#'## Fixed-effect (co)variance matrix
#'vcov(model)
#' 
#' @method coef dosresmeta  
#' @export coef.dosresmeta
#' @S3method coef dosresmeta
#' 
coef.dosresmeta <- function (object, format = c("vector", "matrix"), ...)  {
  coef <- object$coefficients
  format <- match.arg(format, c("vector", "matrix"))
  if (format == "matrix" || is.vector(coef)) 
    return(coef)
  names <- paste(rep(colnames(coef), each = nrow(coef)), rep(rownames(coef), 
                                                             ncol(coef)), sep = "")
  coef <- as.numeric(coef)
  names(coef) <- names
  return(coef)
}

#' @rdname coef.dosresmeta
#' @method vcov dosresmeta
#' @export vcov.dosresmeta
#' @S3method vcov dosresmeta
#' 
vcov.dosresmeta <- function (object, ...) 
{
  return(object$vcov)
}