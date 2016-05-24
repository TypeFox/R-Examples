#' Sample size determination for testing Pearson's correlation coefficient
#'
#' This function performs sample size computation for testing Pearson's correlation coefficient
#' based on precision requirements (i.e., type-I-risk, type-II-risk and an effect size).
#'
#' @param rho            a number indicating the correlation coefficient under the null hypothesis, \eqn{\rho}.0.
#' @param delta          minimum difference to be detected, \eqn{\delta}.
#' @param alternative    a character string specifying the alternative hypothesis,
#'                       must be one of "two.sided" (default), "greater" or "less".
#' @param alpha          type-I-risk, \eqn{\alpha}.
#' @param beta           type-II-risk, \eqn{\beta}.
#' @param output         logical: if \code{TRUE}, output is shown.
#'
#' @author
#' Takuya Yanagida \email{takuya.yanagida@@univie.ac.at},
#'
#' @seealso
#' \code{\link{seqtest.cor}}, \code{\link{size.mean}}, \code{\link{size.prop}}, \code{\link{print.size}}
#'
#' @references
#' Rasch, D., Kubinger, K. D., & Yanagida, T. (2011). \emph{Statistics in psychology - Using R and SPSS}.
#' New York: John Wiley & Sons.
#'
#' Rasch, D., Pilz, J., Verdooren, L. R., & Gebhardt, G. (2011).
#' \emph{Optimal experimental design with R}. Boca Raton: Chapman & Hall/CRC.
#'
#' @return Returns an object of class \code{size} with following entries:
#'
#' \tabular{ll}{
#'   \code{call}      \tab function call \cr
#'   \code{type}      \tab type of the test (i.e., correlation coefficient) \cr
#'   \code{spec}      \tab specification of function arguments \cr
#'   \code{res}       \tab list with the result, i.e., optimal sample size \cr
#' }
#'
#' @export
#'
#' @examples
#'
#' #--------------------------------------
#  # Two-sided test
#' # H0: rho = 0.3, H1: rho != 0.3
#' # alpha = 0.05, beta = 0.2, delta = 0.2
#'
#' size.cor(rho = 0.3, delta = 0.2, alpha = 0.05, beta = 0.2)
#'
#' #--------------------------------------
#  # One-sided test
#' # H0: rho <= 0.3, H1: rho > 0.3
#' # alpha = 0.05, beta = 0.2, delta = 0.2
#'
#' size.cor(rho = 0.3, delta = 0.2, alternative = "greater", alpha = 0.05, beta = 0.2)
size.cor <- function(rho = NULL, delta,
                     alternative = c("two.sided", "less", "greater"),
                     alpha = 0.05, beta = 0.1, output = TRUE) {

  #-----------------------------------------------------------------------------------
  # Input check

  if (delta <= 0) {

    stop("Argument delta out of bound, specify a value > 0")

  }

  ###

  if (is.null(rho)) {

    rho <- 0

  }

  ###

  if (rho <= -1 || rho >= 1) {

    stop("Argument rho out of bound, specify a value between -1 and 1")

  }

  ###

  if (!all(alternative %in% c("two.sided", "less", "greater"))) {

    stop("Argument alternative should be \"two.sided\", \"less\" or \"greater\"")

  }

  ###

  if (alpha <= 0 || alpha >= 1) {

    stop("Argument alpha out of bound, specify a value between 0 and 1")

  }

  ###

  if (beta <= 0 || beta >= 1) {

    stop("Argument beta out of bound, specify a value between 0 and 1")

  }

  #-----------------------------------------------------------------------------------

  # two- or one-sided test
  alternative <- ifelse(all(c("two.sided", "less", "greater") %in% alternative), "two.sided", alternative)

  if (alternative == "two.sided") {

    if ((rho + delta) >= 1 || (rho - delta) <= -1) {

      stop("Value (rho + delta) or (rho - delta) out of bound")

    }

  } else {

    if (alternative == "less") {

      if ((rho - delta) <= -1) {

        stop("Value (rho - delta) out of bound")

      }

    } else {

      if ((rho + delta) >= 1) {

        stop("Value (rho + delta) out of bound")

      }

    }

  }

  #-----------------------------------------------------------------------------------
  # Main function

  alternative <- ifelse(all(c("two.sided", "less", "greater") %in% alternative), "two.sided", alternative)
  side <- switch(alternative, two.sided = 2, less = 1, greater = 1)

  rho.0 <- rho
  rho.1 <- switch(alternative, two.sided = rho.0 + delta, less = rho.0 - delta, greater = rho.0 + delta)

  n <- 3 + 4 * ((qnorm(1 - alpha / side) + qnorm(1 - beta)) / (log((1 + rho.1) / (1 - rho.1)) - log((1 + rho.0) / (1 - rho.0))))^2

  #-----------------------------------------------------------------------------------
  # Return object

  object <- list(call = match.call(),
                 type = "cor",
                 spec = list(delta = delta, rho = rho, alternative = alternative, alpha = alpha, beta = beta),
                 res = list(n = n))

  class(object) <- "size"

  #-----------------------------------------------------------------------------------
  # Output

  if (output == TRUE) { print(object) }

  return(invisible(object))

}
