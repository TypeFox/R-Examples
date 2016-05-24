#' Sequential triangular test for the arithmetic mean
#'
#' This function performs the sequential triangular test for the arithmetic mean
#' in one- or two-samples
#'
#' For the one-sample test, arguments \code{x}, \code{mu} and the minimum difference to be detected
#' has to be specified (i.e., argument \code{y} must not be specified). For the two-sample test, arguments
#' \code{x}, \code{y}, and the minimum difference to be detected has to be specified. There are two options to specify
#' the minimum difference to be detected: (1) using arguments \code{mu}, \code{sigma} and \code{delta} or
#' (2) using arguments \code{mu} and \code{theta}.
#' Note that it is not a requirement to know sigma in advance, i.e., theta can be specified directly. For example,
#' \code{theta = 1} specifies a relative minimum difference to be detected of one standard deviation.
#'
#' In order to specify a one-sided test, argument \code{alternative} has to be used (i.e., two-sided tests are conducted by default).
#' For the one-sample test, \code{alternative = "less"} specifies the null hypothesis, H0: \eqn{\mu} >= \eqn{\mu}.0
#' and the alternative hypothesis, H1: \eqn{\mu} < \eqn{\mu}.0; \code{alternative = "greater"} specifies the
#' null hypothesis, H0: \eqn{\mu} <= \eqn{\mu}.0 and the alternative hypothesis, H1: \eqn{\mu} > \eqn{\mu}.0.
#' For the two-sample test \code{alternative = "less"} specifies the null hypothesis, H0: \eqn{\mu}.1 >= \eqn{\mu}.2
#' and the alternative hypothesis, H1: \eqn{\mu}.1 < \eqn{\mu}.2; \code{alternative = "greater"} specifies
#' the null hypothesis, H0: \eqn{\mu}.1 <= \eqn{\mu}.2 and the alternative hypothesis, H1: \eqn{\mu}.1 > \eqn{\mu}.2.
#'
#' The main characteristic of the sequential triangular test is that there is no fixed sample size given
#' in advance. That is, for the most recent sampling point, one has to decide whether
#' sampling has to be continued or either the null- or the alternative hypothesis can be
#' accepted given specified precision requirements (i.e. type-I-risk, type-II-risk and a minimum difference to be detected).
#' The (cumulative) test statistic \code{Z.m} on a Cartesian coordinate system produces a "sequential path"
#' on a continuation area as a triangle. As long as the statistic remains within that triangle,
#' additional data have to be sampled. If the path touches or exceeds the borderlines of the triangle,
#' sampling is completed. Depending on the particular borderline, the null-hypothesis is either
#' accepted or rejected.
#'
#' @param x              initial data for group x, at least one entry.
#' @param y              initial data for group y, at least one entry for a two-sample test.
#' @param mu             a number indicating the true value of the mean in case of the one-sample test, \eqn{\mu}.0.
#' @param alternative    a character string specifying the alternative hypothesis,
#'                       must be one of "two.sided" (default), "greater" or "less".
#' @param sigma          standard deviation in the population, \eqn{\sigma}.
#' @param delta          absolute minimum difference to be detected, \eqn{\delta}.
#' @param theta          relative minimum difference to be detected, \eqn{\theta}.
#' @param alpha          type-I-risk, \eqn{\alpha}.
#' @param beta           type-II-risk, \eqn{\beta}.
#' @param output         logical: if \code{TRUE}, output is shown.
#' @param plot           logical: if \code{TRUE}, a plot is generated.
#'
#' @author
#' Takuya Yanagida \email{takuya.yanagida@@univie.ac.at},
#'
#' @seealso
#' \code{\link{update.seqtest}}, \code{\link{seqtest.prop}}, \code{\link{seqtest.cor}},
#' \code{\link{print.seqtest}}, \code{\link{plot.seqtest}}, \code{\link{descript}}
#'
#' @references
#' Rasch, D., Pilz, J., Verdooren, L. R., & Gebhardt, G. (2011).
#' \emph{Optimal experimental design with R}. Boca Raton: Chapman & Hall/CRC.
#'
#' Rasch, D., Kubinger, K. D., & Yanagida, T. (2011). \emph{Statistics in psychology - Using R and SPSS}.
#' New York: John Wiley & Sons.
#'
#' @return
#' Returns an object of class \code{seqtest}, to be used for later update steps. The object has
#' following entries:
#'
#' \tabular{ll}{
#'   \code{call}      \tab function call \cr
#'   \code{type}      \tab type of the test (i.e., arithmetic mean) \cr
#'   \code{spec}      \tab specification of function arguments \cr
#'   \code{tri}       \tab specification of the triangular \cr
#'   \code{dat}       \tab data \cr
#'   \code{res}       \tab list with results \cr
#' }
#' @export
#'
#' @examples
#'
#' #--------------------------------------
#' # Two-sided one-sample test
#' # H0: mu = 50, H1: mu != 50
#' # alpha = 0.05, beta = 0.2, theta = 0.5
#'
#' seq.obj <- seqtest.mean(56, mu = 50, theta = 0.5,
#'                         alpha = 0.05, beta = 0.2, plot = TRUE)
#'
#' # alternative specifiation using sigma and delta
#' seq.obj <- seqtest.mean(56, mu = 50, sigma = 10, delta = 5,
#'                         alpha = 0.05, beta = 0.2, plot = TRUE)
#'
#' seq.obj <- update(seq.obj, x = c(54, 52, 46, 49))
#' seq.obj <- update(seq.obj, x = c(46, 49, 51, 45))
#' seq.obj <- update(seq.obj, x = c(51, 42, 50, 53))
#' seq.obj <- update(seq.obj, x = c(50, 53, 49, 53))
#'
#' #--------------------------------------
#' # One-sided one-sample test
#' # H0: mu <= 50, H1: mu > 50
#' # alpha = 0.05, beta = 0.2, theta = 0.5
#'
#' seq.obj <- seqtest.mean(c(56, 53), mu = 50, alternative = "greater",
#'                         theta = 0.5, alpha = 0.05, beta = 0.2, plot = TRUE)
#'
#' # alternative specifiation using sigma and delta
#' seq.obj <- seqtest.mean(c(56, 53), mu = 50, alternative = "greater",
#'                         sigma = 10, delta = 5, alpha = 0.05, beta = 0.2, plot = TRUE)
#'
#' seq.obj <- update(seq.obj, x = c(67, 52, 48, 59))
#' seq.obj <- update(seq.obj, x = c(53, 57, 54, 62))
#' seq.obj <- update(seq.obj, x = 58)
#'
#' #--------------------------------------
#' # Two-sided two-sample test
#' # H0: mu.1 = mu.2, H1: mu.1 != mu.2
#' # alpha = 0.01, beta = 0.1, theta = 1
#'
#' seq.obj <- seqtest.mean(53, 45, theta = 1,
#'                         alpha = 0.01, beta = 0.01, plot = TRUE)
#'
#' # alternative specifiation using sigma and delta
#' seq.obj <- seqtest.mean(57, 45, sigma = 10, delta = 10,
#'                         alpha = 0.01, beta = 0.01, plot = TRUE)
#'
#' seq.obj <- update(seq.obj, x = c(58, 54, 56), y = c(45, 41, 42))
#' seq.obj <- update(seq.obj, x = c(56, 50, 49), y = c(42, 45, 50))
#' seq.obj <- update(seq.obj, x = c(62, 57, 59))
#' seq.obj <- update(seq.obj, y = c(41, 39, 46))
#' seq.obj <- update(seq.obj, x = 67)
#' seq.obj <- update(seq.obj, y = 40)
#' seq.obj <- update(seq.obj, y = 36)
#'
#' #--------------------------------------
#' # One-sided two-sample test
#' # H0: mu.1 <= mu.2, H1: mu.1 > mu.2
#' # alpha = 0.01, beta = 0.1, theta = 1
#'
#' seq.obj <- seqtest.mean(53, 45, alternative = "greater", theta = 1,
#'                         alpha = 0.01, beta = 0.01, plot = TRUE)
#'
#' # alternative specifiation using sigma and delta
#' seq.obj <- seqtest.mean(57, 45, alternative = "greater",sigma = 10, delta = 10,
#'                         alpha = 0.01, beta = 0.01, plot = TRUE)
#'
#' seq.obj <- update(seq.obj, x = c(58, 54, 56), y = c(45, 41, 42))
#' seq.obj <- update(seq.obj, x = c(56, 50, 49), y = c(42, 45, 50))
#' seq.obj <- update(seq.obj, x = c(62, 57, 59))
#' seq.obj <- update(seq.obj, y = c(41, 39, 46))
seqtest.mean <- function(x, y = NULL, mu = NULL,
                         alternative = c("two.sided", "less", "greater"),
                         sigma = NULL, delta = NULL, theta = NULL, alpha = 0.05, beta = 0.1,
                         output = TRUE, plot = FALSE) {

  #-----------------------------------------------------------------------------------
  # Input Check

  if (is.null(y)) {

    if (is.null(mu)) {

      stop("Specify mu for the population mean under the null hypothesis")

    }

  }

  ###

  if (!is.null(y) && !is.null(mu)) {

    stop("Specify x and mu for a one-sample test or x and y for a two-sample test")

  }


  ###

  if (!all(alternative %in% c("two.sided", "less", "greater"))) {

    stop("Argument alternative should be \"two.sided\", \"less\" or \"greater\"")

  }

  ###

  if (is.null(theta)) {

    if (is.null(sigma) || is.null(delta) ) {

      stop("Specify sigma and delta or theta for the minimum difference to be detected")

    }

  }

  ###

  if (!is.null(sigma) && !is.null(delta) && !is.null(theta)) {

    stop("Specify sigma and delta theta for the minimum difference to be detected")

  }

  ###

  if (!is.null(sigma)) {

    if (sigma <= 0) {

      stop("Argument sigma out of bound, specify a value > 0")

    }

  }

  ###

  if (!is.null(delta)) {

    if (delta <= 0) {

      stop("Argument delta out of bound, specify a value > 0")

    }

  }

  ###

  if (!is.null(theta)) {

    if (theta <= 0) {

      stop("Argument theta out of bound, specify a value > 0")

    }

  }

  ###

  if (alpha <= 0 || alpha >= 1) {

    stop("Argument alpha out of bound, specify a value between 0 and 1")

  }

  ###

  if (beta <= 0 || beta >= 1) {

    stop("Argument beta out of bound, specify a value between 0 and 1")

  }

  ###

  if (!is.null(theta)) {

    if (!is.null(delta)) {

      warning("Specified argument delta ignored because theta was specified")

      delta <- NULL

    }

    if (!is.null(sigma)) {

      warning("Specified argument sigma ignored because theta was specified")

      sigma <- NULL

    }

  }

  #-----------------------------------------------------------------------------------

  # sigma
  variance <- ifelse(!is.null(sigma), "known", "unknown")

  # one or two sample
  sample <- ifelse(is.null(y), "one.sample", "two.sample")

  # two- or one-sided test
  alternative <- ifelse(all(c("two.sided", "less", "greater") %in% alternative), "two.sided", alternative)

  #-----------------------------------------------------------------------------------
  # Main function

  ifelse(alternative == "two.sided", u.1a <- qnorm(1 - alpha / 2), u.1a <- qnorm(1 - alpha))
  u.1b <- qnorm(1 - beta)

  # two-sided
  if (alternative == "two.sided") {

    theta1 <- ifelse(variance == "known", -delta / sigma, -theta)
    theta2 <- ifelse(variance == "known",  delta / sigma,  theta)

    a1 <- (1 + u.1b / u.1a) * log(1 / (2 * alpha / 2)) / theta1
    a2 <- (1 + u.1b / u.1a) * log(1 / (2 * alpha / 2)) / theta2

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

      theta <- ifelse(variance == "known", -delta / sigma, -theta)

    } else {

      theta <- ifelse(variance == "known",  delta / sigma, theta)

    }

    a <- (1 + u.1b / u.1a) * log(1 / (2 * alpha)) / theta
    b <- theta / (2 * (1 + u.1b / u.1a))

    V.max <- a / b
    Z.max <- 2 * a

  }

  #-----------------------------------------------------------------------------------
  # Return object

  # one-sample
  if (sample == "one.sample") {

    if (alternative == "two.sided") {

      object <- list(call = match.call(),
                     type = "mean",
                     spec = list(mu = mu, alternative = alternative, sample = sample,
                                 delta = delta, variance = variance, sigma = sigma, theta = theta,
                                 alpha = alpha, beta = beta),
                     tri = list(a1 = a1, a2 = a2, b1 = b1, b2 = b2,
                                u.1a = u.1a, u.1b = u.1b,
                                V.max = V.max, Z.max = Z.max,
                                intersec = intersec),
                     dat = list(x = NULL, n = NULL),
                     res = list(V.m = NULL, Z.m = NULL, decision = "continue", step = 0))

    } else {

      object <- list(call = match.call(),
                     type = "mean",
                     spec = list(mu = mu, alternative = alternative, sample = sample,
                                 delta = delta, variance = variance, sigma = sigma, theta = theta,
                                 alpha = alpha, beta = beta),
                     tri = list(a = a, b = b,
                                u.1a = u.1a, u.1b = u.1b,
                                V.max = V.max, Z.max = Z.max),
                     dat = list(x = NULL, n = NULL),
                     res = list(V.m = NULL, Z.m = NULL, decision = "continue", step = 0))

    }

  # two-sample
  } else {

    if (alternative == "two.sided") {

      object <- list(call = match.call(),
                     type = "mean",
                     spec = list(mu = mu, alternative = alternative, sample = sample,
                                 delta = delta, variance = variance, sigma = sigma, theta = theta,
                                 alpha = alpha, beta = beta),
                     tri = list(a1 = a1, a2 = a2, b1 = b1, b2 = b2,
                                u.1a = u.1a, u.1b = u.1b,
                                V.max = V.max, Z.max = Z.max,
                                intersec = intersec),
                     dat = list(x = NULL, y = NULL, n.1 = NULL, n.2 = NULL),
                     res = list(V.m = NULL, Z.m = NULL, decision = "continue", step = 0))

    } else {

      object <- list(call = match.call(),
                     type = "mean",
                     spec = list(mu = mu, alternative = alternative, sample = sample,
                                 delta = delta, variance = variance, sigma = sigma, theta = theta,
                                 alpha = alpha, beta = beta),
                     tri = list(a = a, b = b,
                                u.1a = u.1a, u.1b = u.1b,
                                V.max = V.max, Z.max = Z.max),
                     dat = list(x = NULL, y = NULL, n.1 = NULL, n.2 = NULL),
                     res = list(V.m = NULL, Z.m = NULL, decision = "continue", step = 0))

    }

  }

  class(object) <- "seqtest"

  #-----------------------------------------------------------------------------------
  # Initial data

  print.step <- 0

  # one-sample
  if (object$spec$sample == "one.sample") {

    print.max <- length(x)

    for (x.i in x) {

      print.step <- print.step + 1

      object <- internal.seqtest.mean(object, x = x.i, initial = TRUE,
                                      print.step = print.step, print.max = print.max, output = output, plot = plot)

      if (object$res$decision != "continue") break

    }

  # two-sample
  } else {

    print.max <- length(c(x, y)) - 1

    print.step <- print.step + 1

    object <- internal.seqtest.mean(object, x = x[1], y = y[1], initial = TRUE,
                                    print.step = print.step, print.max = print.max, output = output, plot = plot)

    x.seq <- seq_along(x)
    y.seq <- seq_along(y)

    # length(x) > 1 | length(y) > 1
    if (max(x.seq, y.seq) > 1) {

      xy.sum <- sum(x.seq %in% y.seq)

      if (xy.sum > 1) {

        for (i in 2:xy.sum) {

          print.step <- print.step + 1

          object <- internal.seqtest.mean(object, x = x[i], initial = TRUE,
                                          print.step = print.step, print.max = print.max, output = output, plot = plot)

          if (object$res$decision != "continue") break

          ###

          print.step <- print.step + 1

          object <- internal.seqtest.mean(object, y = y[i], initial = TRUE,
                                          print.step = print.step, print.max = print.max, output = output, plot = plot)

          if (object$res$decision != "continue") break

        }

      }

      # length(x) > length(y)
      if (length(x) > length(y) && object$res$decision == "continue") {

        for (i in (xy.sum + 1):length(x)) {

          print.step <- print.step + 1

          object <- internal.seqtest.mean(object, x = x[i], initial = TRUE,
                                          print.step = print.step, print.max = print.max, output = output, plot = plot)

          if (object$res$decision != "continue") break

        }

      }

      # length(x) < length(y)
      if (length(x) < length(y) && object$res$decision == "continue") {

        for (i in (xy.sum + 1):length(y)) {

          print.step <- print.step + 1

          object <- internal.seqtest.mean(object, y = y[i], initial = TRUE,
                                          print.step = print.step, print.max = print.max, output = output, plot = plot)

          if (object$res$decision != "continue") break

        }

      }

    }

  }

  #-----------------------------------------------------------------------------------

  return(invisible(object))

}
