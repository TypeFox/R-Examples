#' @title Printing a PPSfit Object
#' @description It prints its argument (typically from \code{PPS.fit()}), returning some of the most important aspects.
#'
#' @param x a \code{PPSfit} object.
#' @param digits the number of digits to be printed.
#' @param \dots other arguments.
#' 
#' @references Sarabia, J.M and Prieto, F. (2009). The Pareto-positive stable distribution: A new descriptive model for city size data, \emph{Physica A: Statistical Mechanics and its Applications}, \bold{388}(19), 4179-4191.
#' 
#' @seealso \code{\link{PPS.fit}}
#'
#' @examples
#' x <- rPPS(50, 1.2, 100, 2.3)
#' fit <- PPS.fit(x)
#' print(fit)

#' @export
print.PPSfit <-
  function (x, digits = max(3, getOption("digits") - 3), ...) 
  {
    if (!class(x) == "PPSfit") {
      stop("Object must belong to class PPS")
    }
    cat("\nData:     ", x$obsName, "\n")
    cat("\n")
    if (!is.null(x$sigma)){
      cat("Sigma:\n")
      print.default(format(x$sigma, digits = digits), print.gap = 2, quote = FALSE)
    }
    cat("\n")
    cat("Parameter estimates:\n")
    print.default(format(x$estimate, digits = digits), print.gap = 2, quote = FALSE)
    cat("\n")
    if (is.null(x$Pareto)) x$Pareto <- FALSE
    if (x$Pareto == TRUE){
      cat("Pareto:\n")
      print.default(TRUE)
    }
    cat("\n")
    cat("Log-likelihood:\n")
    print.default(format(x$loglik, digits = digits), print.gap = 2, quote = FALSE)
    cat("\n")
    cat("Sample size:\n")
    print.default(format(x$n, digits = digits), print.gap = 2, quote = FALSE)
    cat("\n")
    cat("\nEstimation method:     ", x$estim.method, "\n") 
    invisible(x)
  }