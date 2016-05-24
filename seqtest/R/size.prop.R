#' Sample size determination for testing the proportion
#'
#' This function performs sample size computation for the one-sample and two-sample test for proportions
#' based on precision requirements (i.e., type-I-risk, type-II-risk and an effect size).
#'
#' @param pi             a number indicating the true value of the probability under the null hypothesis (one-sample test), \eqn{\pi}.0
#'                       or a number indicating the true value of the probability in group 1 (two-sample test), \eqn{\pi}.1.
#' @param delta          minimum difference to be detected, \eqn{\delta}.
#' @param sample         a character string specifying one- or two-sample proportion test,
#'                       must be one of "two.sample" (default) or "one.sample".
#' @param alternative    a character string specifying the alternative hypothesis,
#'                       must be one of "two.sided" (default), "less" or "greater".
#' @param alpha          type-I-risk, \eqn{\alpha}.
#' @param beta           type-II-risk, \eqn{\beta}.
#' @param correct        a logical indicating whether continuity correction should be applied.
#' @param output         logical: if \code{TRUE}, output is shown.
#'
#' @author
#' Takuya Yanagida \email{takuya.yanagida@@univie.ac.at},
#'
#' @seealso
#' \code{\link{seqtest.prop}}, \code{\link{size.mean}}, \code{\link{size.cor}}, \code{\link{print.size}}
#'
#' @references
#' Fleiss, J. L., Levin, B., & Paik, M. C. (2003). \emph{Statistical methods for rates and proportions} (3rd ed.).
#' New York: John Wiley & Sons.
#'
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
#'   \code{type}      \tab type of the test (i.e., proportion) \cr
#'   \code{spec}      \tab specification of function arguments \cr
#'   \code{res}       \tab list with the result, i.e., optimal sample size \cr
#' }
#'
#' @export
#'
#' @examples
#'
#' #--------------------------------------
#' # Two-sided one-sample test
#' # H0: pi = 0.5, H1: pi != 0.5
#' # alpha = 0.05, beta = 0.2, delta = 0.2
#'
#' size.prop(pi = 0.5, delta = 0.2, sample = "one.sample",
#'           alternative = "two.sided", alpha = 0.05, beta = 0.2)
#'
#' #--------------------------------------
#' # One-sided one-sample test
#' # H0: pi <= 0.5, H1: pi > 0.5
#' # alpha = 0.05, beta = 0.2, delta = 0.2
#'
#' size.prop(pi = 0.5, delta = 0.2, sample = "one.sample",
#'           alternative = "less", alpha = 0.05, beta = 0.2)
#'
#' #--------------------------------------
#' # Two-sided two-sample test
#' # H0: pi.1 = pi.2 = 0.5, H1: pi.1 != pi.2
#' # alpha = 0.01, beta = 0.1, delta = 0.2
#'
#' size.prop(pi = 0.5, delta = 0.2, sample = "two.sample",
#'           alternative = "two.sided", alpha = 0.01, beta = 0.1)
#'
#' #--------------------------------------
#' # One-sided two-sample test
#' # H0: pi.1 <=  pi.1 = 0.5, H1: pi.1 > pi.2
#' # alpha = 0.01, beta = 0.1, delta = 0.2
#'
#' size.prop(pi = 0.5, delta = 0.2, sample = "two.sample",
#'           alternative = "greater", alpha = 0.01, beta = 0.1)
size.prop <- function(pi = NULL, delta, sample = c("two.sample", "one.sample"),
                      alternative = c("two.sided", "less", "greater"),
                      alpha = 0.05, beta = 0.1, correct = FALSE, output = TRUE) {

  #-----------------------------------------------------------------------------------
  # Input check

  if (delta <= 0) {

    stop("Argument theta out of bound, specify a value > 0")

  }

  ###

  if (is.null(pi)) {

    pi <- 0.5

  }

  ###

  if (pi >= 1 || pi <= 0) {

    stop("Argument pi out of bound, specify a value between 0 and 1")

  }

  ###

  if (!all(sample %in% c("two.sample", "one.sample"))) {

    stop("Argument sample should be \"two.siample\" or \"one.sample\"")

  }


  ###

  if (!all(alternative %in% c("two.sided", "less", "greater"))) {

    stop("Argument alternative should be \"two.sided\", \"less\", or \"greater\"")

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

  # one or two sample
  sample <- ifelse(all(c("two.sample", "one.sample") %in% alternative), "two.sample", sample)

  # two- or one-sided test
  alternative <- ifelse(all(c("two.sided", "less", "greater") %in% alternative), "two.sided", alternative)

  ###

  if (alternative == "two.sided") {

    if ((pi + delta) >= 1 || (pi - delta) <= 0) {

      stop("Value (pi + delta) or (pi - delta) out of bound")

    }

  } else {

    # one-sample
    if (sample == "one.sample") {

      if (alternative == "less") {

        if ((pi - delta) <= 0) {

          stop("Value (pi - delta) out of bound")

        }

      } else {


        if ((pi + delta) >= 1) {

          stop("Value (pi + delta) out of bound")

        }

      }

    # two-sample
    } else {

      if (alternative == "less") {

        if ((pi + delta) >= 1) {

          stop("Value (pi + delta) out of bound")

        }

      } else {

        if ((pi - delta) <= 0) {

          stop("Value (pi - delta) out of bound")

        }

      }

    }

  }

  #-----------------------------------------------------------------------------------
  # Main function

  side <- switch(alternative, two.sided = 2, less = 1, greater = 1)

  # two-sample
  if (sample == "two.sample") {

    pi.1 <- pi
    pi.2 <- switch(alternative, two.sided = pi.1 + delta, less = pi.1 + delta, greater = pi.1 - delta)

    p.body <- quote({
      pnorm((sqrt(n) * abs(pi.1 - pi.2) - (qnorm(1 - alpha / side) * sqrt((pi.1 + pi.2) * (1 - (pi.1 + pi.2) / 2)))) / sqrt(pi.1 * (1 - pi.1) + pi.2 * (1 - pi.2)))
    })

    n <- uniroot(function(n) eval(p.body) - (1 - beta), c(1, 1e+07))$root

    if (correct == TRUE) {

      n <- ceiling(n)
      n <- (n / 4) * (1 + sqrt(1 + 4 / (n * delta)))^2

    }

  # one-sample
  } else {

    pi.0 <- pi
    pi.1 <- switch(alternative, two.sided = pi.0 + delta, less = pi.0 - delta, greater = pi.0 + delta)

    n <- ((qnorm(1 - alpha / side) * sqrt(pi.0 * (1 - pi.0)) + qnorm(1 - beta) * sqrt(pi.1 * (1 - pi.1))) / (pi.1 - pi.0))^2

    if (correct == TRUE) {

      n <- ceiling(n)
      n <- n + 1 / (qnorm(1 - alpha / side) * sqrt(pi.0 * (1 - pi.0) / n) + qnorm(1 - beta) * sqrt(pi.1 * (1 - pi.1) / n))

    }

  }

  #-----------------------------------------------------------------------------------
  # Return object

  object <- list(call = match.call(),
                 type = "prop",
                 spec = list(delta = delta, pi = pi, sample = sample, alternative = alternative,
                             alpha = alpha, beta = beta, correct = correct),
                 res = list(n = n))

  class(object) <- "size"

  #-----------------------------------------------------------------------------------
  # Output

  if (output == TRUE) { print(object) }

  return(invisible(object))

}
