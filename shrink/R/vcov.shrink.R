#' Calculate Variance-Covariance Matrix of Shrinkage Factors for Objects of Class \code{shrink}
#'
#' This class of objects is returned by the \code{shrink} function. Objects of this
#' class have methods for the functions \code{coef}, \code{predict}, \code{print},
#' \code{summary}, and \code{vcov}.
#'
#' @param object object of class \code{shrink}.
#' @param digits integer, used for number formatting with \code{\link[base]{signif}}().
#' @param ... further arguments.
#'
#' @export
#'
#' @return A matrix of the estimated covariances between the obtained shrinkage factors.
#'
#' @seealso \code{\link{shrink}},\code{\link{coef.shrink}}, \code{\link{predict.shrink}},
#'     \code{\link{print.shrink}}, \code{\link{summary.shrink}}
vcov.shrink <-
  function(object, digits = 6, ...)
{
  if (!inherits(object, "shrink")) stop("'object' is not of class shrink")

  if (object$type == "all") {
    result <- vector(mode = "list", length = 2)
    names(result) <- c("global", "parameterwise")

    cat("GLOBAL SHRINKAGE\n")
    result[[1]] <- print(signif (object$global$ShrinkageFactorsVCOV, digits = digits))

    cat("\n\nPARAMETERWISE SHRINKAGE\n")
    result[[2]] <- print(signif (object$parameterwise$ShrinkageFactorsVCOV,
                                 digits = digits))

    if (!is.null(object$join)) {
      cat("\n\nPARAMETERWISE SHRINKAGE WITH JOIN OPTION\n")
      result$joint <- print(signif (object$joint$ShrinkageFactorsVCOV, digits = digits))
      if (!is.null(object$join)) cat("\njoint shrinkage was requested for:",
                                     sapply(object$join, function(x) paste(x, collapse="+")))
    }
  } else
    result <- print(signif (object$ShrinkageFactorsVCOV, digits = digits))

  invisible(result)
}
