#' Print Method for Objects of Class \code{shrink}
#'
#' This class of objects is returned by the \code{shrink} function. Objects of this
#' class have methods for the functions \code{coef}, \code{predict}, \code{print},
#' \code{summary}, and \code{vcov}.
#'
#' @param x object of class \code{shrink}.
#' @param ... further arguments.
#'
#' @export
#'
#' @seealso \code{\link{shrink}}, \code{\link{coef.shrink}}, \code{\link{predict.shrink}},
#'     \code{\link{summary.shrink}}, \code{\link{vcov.shrink}}
print.shrink <- function(x, ...)
{
  if (!inherits(x, "shrink")) stop("'object' is not of class shrink")

  if (x$type == "all") {
    cat("GLOBAL SHRINKAGE\n", as.vector(x$global$ShrinkageFactors))
    cat("\nShrunken Regression Coefficients:\n")
    print(x$global$ShrunkenRegCoef)

    cat("\n\nPARAMETERWISE SHRINKAGE\n")
    print(x$parameterwise$ShrinkageFactors)
    cat("\nShrunken Regression Coefficients:\n")
    print(x$parameterwise$ShrunkenRegCoef)

    if (!is.null(x$join)) {
      cat("\n\nPARAMETERWISE SHRINKAGE WITH JOIN OPTION\n")
      print(x$joint$ShrinkageFactors)
      cat("\nShrunken Regression Coefficients:\n")
      print(x$joint$ShrunkenRegCoef)
    }
  } else {
    if (is.null(x$join)) {
      cat(paste("Shrinkage Factors (type=", x$type, ", method=", x$method, "):\n",
                sep=""))
    } else {
      cat(paste("Shrinkage Factors (type=", x$type, " with join option, method=",
                x$method, "):\n", sep=""))
    }
    if (x$type %in% "global") { print(as.vector(x$ShrinkageFactors)) } else {
      print(x$ShrinkageFactors) }

    #  if (!is.null(x$postfit)) {
    #    cat("\nPostfit:\n")
    #    print(x$postfit, ...)
    #  } else {
    cat("\nShrunken Regression Coefficients:\n")
    print(x$ShrunkenRegCoef, ...)
    #  }
  }
}
