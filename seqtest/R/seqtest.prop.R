#' Sequential triangular test for the proportion
#'
#' This function performs the sequential triangular test for the proportion
#' in one- or two-samples
#'
#' For the one-sample test, arguments \code{x}, \code{pi}, and \code{delta} has to be specified
#' (i.e., argument \code{y} must not be specified).
#' For the two-sample test, arguments \code{x}, \code{y}, \code{pi}, and \code{delta} has to be specified
#'
#' In order to specify a one-sided test, argument \code{alternative} has to be used (i.e., two-sided tests are conducted by default).
#' For the one-sample test, \code{alternative = "less"} specifies the null hypothesis, H0: \eqn{\pi} >= \eqn{\pi}.0
#' and the alternative hypothesis, H1: \eqn{\pi} < \eqn{\pi}.0; \code{alternative = "greater"} specifies the
#' null hypothesis, H0: \eqn{\pi} <= \eqn{\pi}.0 and the alternative hypothesis, H1: \eqn{\pi} > \eqn{\pi}.0.
#' For the two-sample test \code{alternative = "less"} specifies the null hypothesis, H0: \eqn{\pi}.1 >= \eqn{\pi}.2
#' and the alternative hypothesis, H1: \eqn{\pi}.1 < \eqn{\pi}.2; \code{alternative = "greater"} specifies
#' the null hypothesis, H0: \eqn{\pi}.1 <= \eqn{\pi}.2 and the alternative hypothesis, H1: \eqn{\pi}.1 > \eqn{\pi}.2.
#'
#' The main characteristic of the sequential triangular test is that there is no fixed sample size given
#' in advance. That is, for the most recent sampling point, one has to decide whether
#' sampling has to be continued or either the null- or the alternative hypothesis can be
#' accepted given specified precision requirements (i.e. type-I-risk, type-II-risk and an effect size).
#' The (cumulative) test statistic \code{Z.m} on a Cartesian coordinate system produces a "sequential path"
#' on a continuation area as a triangle. As long as the statistic remains within that triangle,
#' additional data have to be sampled. If the path touches or exceeds the borderlines of the triangle,
#' sampling is completed. Depending on the particular borderline, the null-hypothesis is either
#' accepted or rejected.
#'
#' @param x              initial data for group x, at least one entry.
#' @param y              initial data for group y, at least one entry for a two-sample test.
#' @param pi             a number indicating the true value of the probability of success in group x, \eqn{\pi}.0.
#' @param alternative    a character string specifying the alternative hypothesis,
#'                       must be one of "two.sided" (default), "less" or "greater".
#' @param delta          minimum difference to be detected, \eqn{\delta}.
#' @param alpha          type-I-risk, \eqn{\alpha}.
#' @param beta           type-II-risk, \eqn{\beta}.
#' @param output         logical: if \code{TRUE}, output is shown.
#' @param plot           logical: if \code{TRUE}, a plot is generated.
#'
#' @author
#' Takuya Yanagida \email{takuya.yanagida@@univie.ac.at},
#'
#' @seealso
#' \code{\link{update.seqtest}}, \code{\link{seqtest.mean}}, \code{\link{seqtest.cor}},
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
#'   \code{type}      \tab type of the test (i.e., proportion) \cr
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
#' # H0: pi = 0.5, H1: pi != 0.5
#' # alpha = 0.05, beta = 0.2, delta = 0.2
#'
#' seq.obj <- seqtest.prop(c(1, 1, 0, 1), pi = 0.5, delta = 0.2,
#'                         alpha = 0.05, beta = 0.2, plot = TRUE)
#'
#' seq.obj <- update(seq.obj, x = c(1, 1, 1, 1, 1, 0, 1, 1, 1))
#' seq.obj <- update(seq.obj, x = c(0, 1, 1, 1))
#' seq.obj <- update(seq.obj, x = c(1, 1))
#'
#' #--------------------------------------
#' # One-sided one-sample test
#' # H0: pi <= 0.5, H1: pi > 0.5
#' # alpha = 0.05, beta = 0.2, delta = 0.2
#'
#' seq.obj <- seqtest.prop(c(1, 1, 0, 1), pi = 0.5,
#'                         alternative = "greater", delta = 0.2,
#'                         alpha = 0.05, beta = 0.2, plot = TRUE)
#'
#' seq.obj <- update(seq.obj, x = c(1, 1, 1, 1, 1, 0, 1, 1, 1))
#' seq.obj <- update(seq.obj, x = c(0, 1, 1, 1))
#'
#' #--------------------------------------
#' # Two-sided two-sample test
#' # H0: pi.1 = pi.2 = 0.5, H1: pi.1 != pi.2
#' # alpha = 0.01, beta = 0.1, delta = 0.2
#'
#' seq.obj <- seqtest.prop(1, 0, pi = 0.5, delta = 0.2,
#'                         alpha = 0.01, beta = 0.1, plot = TRUE)
#'
#' seq.obj <- update(seq.obj, x = c(1, 1, 1, 0), y = c(0, 0, 1, 0))
#' seq.obj <- update(seq.obj, x = c(0, 1, 1, 1), y = c(0, 0, 0, 0))
#' seq.obj <- update(seq.obj, x = c(1, 0, 1, 1), y = c(0, 0, 0, 1))
#' seq.obj <- update(seq.obj, x = c(1, 1, 1, 1), y = c(0, 0, 0, 0))
#' seq.obj <- update(seq.obj, x = c(0, 1, 0, 1))
#' seq.obj <- update(seq.obj, y = c(0, 0, 0, 1))
#' seq.obj <- update(seq.obj, x = c(1, 1, 1, 1))
#'
#' #--------------------------------------
#' # One-sided two-sample test
#' # H0: pi.1 <=  pi.1 = 0.5, H1: pi.1 > pi.2
#' # alpha = 0.01, beta = 0.1, delta = 0.2
#'
#' seq.obj <- seqtest.prop(1, 0, pi = 0.5, delta = 0.2,
#'                         alternative = "greater",
#'                         alpha = 0.01, beta = 0.1, plot = TRUE)
#'
#' seq.obj <- update(seq.obj, x = c(1, 1, 1, 0), y = c(0, 0, 1, 0))
#' seq.obj <- update(seq.obj, x = c(0, 1, 1, 1), y = c(0, 0, 0, 0))
#' seq.obj <- update(seq.obj, x = c(1, 0, 1, 1), y = c(0, 0, 0, 1))
#' seq.obj <- update(seq.obj, x = c(1, 1, 1), y = c(0, 0))
seqtest.prop <- function(x, y = NULL, pi = NULL,
                         alternative = c("two.sided", "less", "greater"),
                         delta, alpha = 0.05, beta = 0.1,
                         output = TRUE, plot = FALSE) {

  #-----------------------------------------------------------------------------------
  # Input Check

  if (any(!x %in% c(0, 1)) || any(!y %in% c(0, 1))) {

    stop("Only 0 and 1 are allowed for x and y")

  }

  ###

  if (is.null(pi)) {

    pi <- 0.5

  }

  ###

  if (pi <= 0 || pi >= 1) {

    stop("Argument pi out of bound, specify a value between 0 and 1")

  }

  ###

  if (!all(alternative %in% c("two.sided", "less", "greater"))) {

    stop("Argument alternative should be \"two.sided\", \"less\" or \"greater\"")

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

  # one or two sample
  sample <- ifelse(is.null(y), "one.sample", "two.sample")

  # two- or one-sided test
  alternative <- ifelse(all(c("two.sided", "less", "greater") %in% alternative), "two.sided", alternative)

  ###

  if (alternative == "two.sided") {

    if ((pi + delta) >= 1 || (pi - delta) <= 0) {

      stop("Value pi + delta or pi - delta out of bound")

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

  ifelse(alternative == "two.sided", u.1a <- qnorm(1 - alpha / 2), u.1a <- qnorm(1 - alpha))
  u.1b <- qnorm(1 - beta)

  # two-sided
  if (alternative == "two.sided") {

    if (sample == "one.sample") {

      theta1 <-  log(((pi - delta) * (1 - pi)) / (pi * (1 - (pi - delta))))
      theta2 <-  log(((pi + delta) * (1 - pi)) / (pi * (1 - (pi + delta))))

    } else {

     theta1 <- -log(((pi + delta) * (1 - pi)) / (pi * (1 - (pi + delta))))
     theta2 <- -log(((pi - delta) * (1 - pi)) / (pi * (1 - (pi - delta))))

    }

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

    if (sample == "one.sample") {

      if (alternative == "less") {

        theta <- log(((pi - delta) * (1 - pi)) / (pi * (1 - (pi - delta))))

      } else {

        theta <- log(((pi + delta) * (1 - pi)) / (pi * (1 - (pi + delta))))

      }

    } else {

      if (alternative == "less") {

        theta <- -log(((pi + delta) * (1 - pi)) / (pi * (1 - (pi + delta))))

      } else {

        theta <- -log(((pi - delta) * (1 - pi)) / (pi * (1 - (pi - delta))))

      }

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
                     type = "prop",
                     spec = list(pi = pi, alternative = alternative, sample = sample,
                                 delta = delta, theta1 = theta1, theta2 = theta2,
                                 alpha = alpha, beta = beta),
                     tri = list(a1 = a1, a2 = a2, b1 = b1, b2 = b2,
                                u.1a = u.1a, u.1b = u.1b,
                                V.max = V.max, Z.max = Z.max,
                                intersec = intersec),
                     dat = list(x = NULL, n = NULL),
                     res = list(V.m = NULL, Z.m = NULL, decision = "continue", step = 0))

    } else {

      object <- list(call = match.call(),
                     type = "prop",
                     spec = list(pi = pi, alternative = alternative, sample = sample,
                                 delta = delta, theta = theta,
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
                     type = "prop",
                     spec = list(pi = pi, alternative = alternative, sample = sample,
                                 delta = delta, theta1 = theta1, theta2 = theta2,
                                 alpha = alpha, beta = beta),
                     tri = list(a1 = a1, a2 = a2, b1 = b1, b2 = b2,
                                u.1a = u.1a, u.1b = u.1b,
                                V.max = V.max, Z.max = Z.max,
                                intersec = intersec),
                     dat = list(x = NULL, y = NULL, n.1 = NULL, n.2 = NULL),
                     res = list(V.m = NULL, Z.m = NULL, decision = "continue", step = 0))

    } else {

      object <- list(call = match.call(),
                     type = "prop",
                     spec = list(pi = pi, alternative = alternative, sample = sample,
                                 delta = delta, theta = theta,
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

      object <- internal.seqtest.prop(object, x = x.i, initial = TRUE,
                                      print.step = print.step, print.max = print.max, output = output, plot = plot)

      if (object$res$decision != "continue") break

    }

  # two-sample
  } else {

    print.max <- length(c(x, y)) - 1

    print.step <- print.step + 1

    object <- internal.seqtest.prop(object, x = x[1], y = y[1], initial = TRUE,
                                    print.step = print.step, print.max = print.max, output = output, plot = plot)

    x.seq <- seq_along(x)
    y.seq <- seq_along(y)

    # length(x) > 1 | length(y) > 1
    if (max(x.seq, y.seq) > 1) {

      xy.sum <- sum(x.seq %in% y.seq)

      if (xy.sum > 1) {

        for (i in 2:xy.sum) {

          print.step <- print.step + 1

          object <- internal.seqtest.prop(object, x = x[i], initial = TRUE,
                                          print.step = print.step, print.max = print.max, output = output, plot = plot)

          if (object$res$decision != "continue") break

          ###

          print.step <- print.step + 1

          object <- internal.seqtest.prop(object, y = y[i], initial = TRUE,
                                          print.step = print.step, print.max = print.max, output = output, plot = plot)

          if (object$res$decision != "continue") break

        }

      }

      # length(x) > length(y)
      if (length(x) > length(y) & object$res$decision == "continue") {

        for (i in (xy.sum + 1):length(x)) {

          print.step <- print.step + 1

          object <- internal.seqtest.prop(object, x = x[i], initial = TRUE,
                                          print.step = print.step, print.max = print.max, output = output, plot = plot)

          if (object$res$decision != "continue") break

        }

      }

      # length(x) < length(y)
      if (length(x) < length(y) & object$res$decision == "continue") {

        for (i in (xy.sum + 1):length(y)) {

          print.step <- print.step + 1

          object <- internal.seqtest.prop(object, y = y[i], initial = TRUE,
                                          print.step = print.step, print.max = print.max, output = output, plot = plot)

          if (object$res$decision != "continue") break

        }

      }

    }

  }

  #-----------------------------------------------------------------------------------

  return(invisible(object))

}
