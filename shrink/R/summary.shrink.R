#' Summary Method for Objects of Class \code{shrink}
#'
#' This class of objects is returned by the \code{shrink} function. Objects of this
#' class have methods for the functions \code{coef}, \code{predict}, \code{print},
#' \code{summary}, and \code{vcov}.
#'
#' @param object an object of class \code{shrink}.
#' @param digits integer, used for number formatting with \code{\link[base]{signif}}().
#' @param ... further arguments.
#'
#' @export
#'
#' @return A matrix with regression coefficients of the orginial fit, corresponding
#'     shrinkage factors and shrunken regression coefficients.
#' @seealso \code{\link{shrink}}, \code{\link{coef.shrink}}, \code{\link{print.shrink}},
#'     \code{\link{predict.shrink}}, \code{\link{vcov.shrink}}
summary.shrink <-
  function(object, digits = 6, ...)
{
  if (!inherits(object, "shrink")) stop("'object' is not of class shrink")

  cat("Call:\n")
  print(object$fit$call)
  print(object$call)

  cat("\nCoefficients:\n")
  if (object$type == "all") {
    #    if (inherits(object$fit, "coxph")) interceptG <- interceptP <- NULL else
    if (attr(terms(object$fit), "intercept") && !inherits(object$fit, "coxph")) {
      interceptG <- object$global$ShrunkenRegCoef[1L] / object$fit$coefficients[1L]
      interceptP <- object$parameterwise$ShrunkenRegCoef[1L] / object$fit$coefficients[1L]
      nb <- length(object$fit$coefficients) - 1L
    } else {
      interceptP <- interceptG <- NULL
      nb <- length(object$fit$coefficients)
    }

    output <- matrix(cbind(if (inherits(object$fit, "mfp")) object$fit$fit$coefficients else object$fit$coefficients,
                           c(interceptG, rep(object$global$ShrinkageFactors, times = nb)),
                           object$global$ShrunkenRegCoef,
                           c(interceptP, object$parameterwise$ShrinkageFactors),
                           object$parameterwise$ShrunkenRegCoef),
                     ncol = 5, dimnames=list(names(object$global$ShrunkenRegCoef),
                                             c("Estimate", "gShrinkageFactor", "gShrunkenEst.",
                                               "pShrinkageFactor", "pShrunkenEst.")))

    if (!is.null(object$join)) {
      if (!inherits(object$fit, "coxph") && attr(terms(object$fit), "intercept")) interceptJ <- object$joint$ShrunkenRegCoef[1L] / object$fit$coefficients[1L] else interceptJ <- NULL
        output <- cbind(output,
                        'jShrinkageFactor' = c(interceptJ, object$joint$ShrinkageFactors),
                        'j ShrunkenEst.' = object$joint$ShrunkenRegCoef)
    }
  } else {
    if (object$type == "global") {
      if (inherits(object$fit, "coxph")) {
        interceptG <- NULL
        nb <- length(object$fit$coefficients)
      } else
        if (attr(terms(object$fit), "intercept")) {
          interceptG <- object$ShrunkenRegCoef[1L] / object$fit$coefficients[1L]
          nb <- length(object$fit$coefficients) - 1L
        } else {
          interceptG <- NULL
          nb <- length(object$fit$coefficients)
        }

      output <- matrix(cbind(if (inherits(object$fit, "mfp")) object$fit$fit$coefficients else object$fit$coefficients,
                             c(interceptG, rep(object$ShrinkageFactors, times = nb)),
                             object$ShrunkenRegCoef),
                       ncol = 3, dimnames=list(names(object$ShrunkenRegCoef),
                                               c("Estimate", paste(substr(object$type, 1L, 1L),
                                                                   c("ShrinkageFactor", "ShrunkenEstimate"), sep = ""))))
    } else {
      if (inherits(object$fit, "coxph")) interceptP <- NULL else
        if (attr(terms(object$fit), "intercept")) interceptP <- object$ShrunkenRegCoef[1L] / object$fit$coefficients[1L] else interceptP <- NULL

        output <- matrix(cbind(if (inherits(object$fit, "mfp")) object$fit$fit$coefficients else object$fit$coefficients,
                               c(interceptP, object$ShrinkageFactors),
                               object$ShrunkenRegCoef),
                         ncol = 3, dimnames=list(names(object$ShrunkenRegCoef),
                                                 c("Estimate", paste(substr(object$type, 1L, 1L),
                                                                     c("ShrinkageFactor", "ShrunkenEstimate"), sep = ""))))
    }
  }

  print(signif (output, digits))

  if (!is.null(object$join)) cat("\njoint shrinkage was requested for:",
                                 sapply(object$join, function(object) paste(object, collapse="+")))
  cat("\n\ng = global shrinkage; p = parameterwise shrinkage; j = parameterwise shrinkage with the joint option")

  invisible(signif (output, digits))
}
