#' Sequential triangular test for Pearson's correlation coefficient
#'
#' This function performs the sequential triangular test for Pearson's correlation coefficient
#'
#' Null and alternative hypothesis is specified using arguments \code{rho} and \code{delta}.
#' Note that the argument k (i.e., number of observations in each sub-sample) has to be specified. At least k = 4 is needed.
#' The optimal value of k should be determined based on statistical simulation using \code{\link{sim.seqtest.cor}} function.
#'
#' In order to specify a one-sided test, argument \code{alternative} has to be used (i.e., two-sided tests are conducted by default).
#' That is, \code{alternative = "less"} specifies the null hypothesis, H0: \eqn{\rho} >= \eqn{\rho}.0 and
#' the alternative hypothesis, H1: \eqn{\rho} < \eqn{\rho}.0; \code{alternative = "greater"} specifies the
#' null hypothesis, H0: \eqn{\rho} <= \eqn{\rho}.0 and the alternative hypothesis, H1: \eqn{\rho} > \eqn{\rho}.0.
#'
#' The main characteristic of the sequential triangular test is that there is no fixed sample size given
#' in advance. That is, for the most recent sampling point, one has to decide whether
#' sampling has to be continued or either the null- or the alternative hypothesis can be
#' accepted given specified precision requirements (i.e. type-I-risk, type-II-risk and an effect size).
#' The sequence of data  pairs must we split into sub-samples of length k >= 4 each.
#' The (cumulative) test statistic \code{Z.m} on a Cartesian coordinate system produces a "sequential path" on a
#' continuation area as a triangle. As long as the statistic remains within that triangle,
#' additional data have to be sampled. If the path touches or exceeds the borderlines of the triangle,
#' sampling is completed. Depending on the particular borderline, the null-hypothesis is either
#' accepted or rejected.
#'
#' @param x              initial data, i.e., Pearson's correlation coefficient in a sub-sample of k observations.
#' @param k              number of observations in each sub-sample.
#' @param rho            a number indicating the correlation coefficient under the null hypothesis, \eqn{\rho}.0.
#' @param alternative    a character string specifying the alternative hypothesis,
#' @param delta          minimum difference to be detected, \eqn{\delta}.
#' @param alpha          type-I-risk, \eqn{\alpha}.
#' @param beta           type-II-risk, \eqn{\beta}.
#' @param output         logical: if \code{TRUE}, output is shown.
#' @param plot           logical: if \code{TRUE}, an initial plot is generated.
#'
#' @author
#' Takuya Yanagida \email{takuya.yanagida@@univie.ac.at},
#'
#' @seealso
#' \code{\link{update.seqtest}}, \code{\link{sim.seqtest.cor}}, \code{\link{seqtest.mean}}, \code{\link{seqtest.prop}},
#' \code{\link{print.seqtest}}, \code{\link{plot.seqtest}}, \code{\link{descript}}
#'
#' @references
#' Schneider, B., Rasch, D., Kubinger, K. D., & Yanagida, T. (2015).
#' A Sequential triangular test of a correlation coefficient's null-hypothesis: 0 \eqn{< \rho \le \rho}0.
#' \emph{Statistical Papers, 56}, 689-699.
#'
#' @return
#' Returns an object of class \code{seqtest}, to be used for later update steps. The object has
#' following entries:
#'
#' \tabular{ll}{
#'   \code{call}      \tab function call \cr
#'   \code{type}      \tab type of the test (i.e., correlation coefficient) \cr
#'   \code{spec}      \tab specification of function arguments \cr
#'   \code{tri}       \tab specification of triangular \cr
#'   \code{dat}       \tab data \cr
#'   \code{res}       \tab list with results \cr
#' }
#' @export
#'
#' @examples
#'
#' #--------------------------------------
#  # Two-sided test
#' # H0: rho = 0.3, H1: rho != 0.3
#' # alpha = 0.05, beta = 0.2, delta = 0.2
#'
#' seq.obj <- seqtest.cor(0.46, k = 14, rho = 0.3, delta = 0.2,
#'                        alpha = 0.05, beta = 0.2, plot = TRUE)
#'
#' seq.obj <- update(seq.obj, c(0.56, 0.76, 0.56, 0.52))
#'
#' #--------------------------------------
#  # One-sided test
#' # H0: rho <= 0.3, H1: rho > 0.3
#' # alpha = 0.05, beta = 0.2, delta = 0.2
#'
#' seq.obj <- seqtest.cor(0.46, k = 14, rho = 0.3,
#'                        alternative = "greater", delta = 0.2,
#'                        alpha = 0.05, beta = 0.2, plot = TRUE)
#'
#' seq.obj <- update(seq.obj, c(0.56, 0.76, 0.66))
seqtest.cor <- function(x, k, rho,
                        alternative = c("two.sided", "less", "greater"),
                        delta, alpha = 0.05, beta = 0.1,
                        output = TRUE, plot = FALSE) {

  #-----------------------------------------------------------------------------------
  # Input Check

  if (any(x < -1) || any(x > 1)) {

    stop("Correlation coefficient out of bound, only values between -1 and 1 are allowed for x")

  }

  ###

  if (k < 4) {

    stop("At least k = 4 data pairs are needed")

  }

  ###

  if (!all(alternative %in% c("two.sided", "less", "greater"))) {

    stop("Argument alternative should be \"two.sided\", \"less\" or \"greater\"")

  }


  ###

  if (rho <= -1 || rho >= 1) {

    stop("Argument rho out of bound, specify a value between -1 and 1")

  }

  ###

  if (delta <= 0) {

    stop("Argument delta out of bound, specify a value > 0")

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

  ifelse(alternative == "two.sided", u.1a <- qnorm(1 - alpha / 2), u.1a <- qnorm(1 - alpha))
  u.1b <- qnorm(1 - beta)

  sd.0 <- sqrt(k - 3) / 2
  z.0 <- log((1 + rho) / (1 - rho)) + rho / (k - 1)

  # two-sided
  if (alternative == "two.sided") {

    z.1 <- log((1 + (rho - delta)) / (1 - (rho - delta))) + (rho - delta) / (k - 1)
    theta1 <- sd.0 * (z.1 - z.0)

    z.1 <- log((1 + (rho + delta)) / (1 - (rho + delta))) + (rho + delta) / (k - 1)
    theta2 <- sd.0 * (z.1 - z.0)

    a1 <- ((1 + u.1b / u.1a) * log(1 / (2 * alpha))) / theta1
    a2 <- ((1 + u.1b / u.1a) * log(1 / (2 * alpha))) / theta2

    b1 <- theta1 / (2 * (1 + u.1b / u.1a))
    b2 <- theta2 / (2 * (1 + u.1b / u.1a))

    V.max <- c(a1 / b1, a2 / b2)
    Z.max <- c(2 * a1, 2 * a2)

    # point of intersection
    f1 <- function(x1) { -a1 + 3 * b1 * x1 }
    f2 <- function(x2) { -a2 + 3 * b2 * x2 }

    intersec <- uniroot(function(z) f1(z) - f2(z), interval = c(0, max(V.max)))$root

  # one-sided
  } else {

    if (alternative == "less") {

      z.1 <- log((1 + (rho - delta)) / (1 - (rho - delta))) + (rho - delta) / (k - 1)

    } else {

      z.1 <- log((1 + (rho + delta)) / (1 - (rho + delta))) + (rho + delta) / (k - 1)

    }

    theta <- sd.0 * (z.1 - z.0)

    a <- ((1 + u.1b / u.1a) * log(1 / (2 * alpha))) / theta
    b <- theta / (2 * (1 + u.1b / u.1a))

    V.max <- a / b
    Z.max <- 2 * a

  }

  #-----------------------------------------------------------------------------------
  # Return object

  if (alternative == "two.sided") {

    object <- list(call = match.call(),
                   type = "cor",
                   spec = list(k = k, rho = rho, alternative = alternative,
                               delta = delta, theta1 = theta1, theta2 = theta2,
                               alpha = alpha, beta = beta),
                   tri = list(a1 = a1, a2 = a2, b1 = b1, b2 = b2,
                              sd.0 = sd.0, z.0 = z.0,
                              u.1a = u.1a, u.1b = u.1b,
                              V.max = V.max, Z.max = Z.max,
                              intersec = intersec),
                   dat = list(x = NULL),
                   res = list(V.m = NULL, Z.m = NULL, decision = "continue", step = 0))

  } else {

    object <- list(call = match.call(),
                   type = "cor",
                   spec = list(k = k, rho = rho, alternative = alternative,
                               delta = delta, theta = theta,
                               alpha = alpha, beta = beta),
                   tri = list(a = a, b = b,
                              sd.0 = sd.0, z.0 = z.0,
                              u.1a = u.1a, u.1b = u.1b,
                              V.max = V.max, Z.max = Z.max),
                   dat = list(x = NULL),
                   res = list(V.m = NULL, Z.m = NULL, decision = "continue", step = 0))

  }

  class(object) <- "seqtest"

  #-----------------------------------------------------------------------------------
  # Initial data

  print.step <- 0
  print.max <- length(x)

  for (x.i in x) {

    print.step <- print.step + 1

    object <- internal.seqtest.cor(object, x = x.i, initial = TRUE,
                                   print.step = print.step, print.max = print.max, output = output, plot = plot)

    if (object$res$decision != "continue") break

  }

  #-----------------------------------------------------------------------------------

  return(invisible(object))

}
