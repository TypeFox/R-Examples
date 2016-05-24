#' Print size object
#'
#' This function prints the \code{size} object
#'
#' @param x           \code{size} object.
#' @param ...         further arguments passed to or from other methods.
#'
#' @author
#' Takuya Yanagida \email{takuya.yanagida@@univie.ac.at}
#'
#' @seealso
#' \code{\link{size.mean}}, \code{\link{size.prop}}, \code{\link{size.cor}}
#'
#' @references
#' Rasch, D., Kubinger, K. D., & Yanagida, T. (2011). \emph{Statistics in psychology - Using R and SPSS}.
#' New York: John Wiley & Sons.
#'
#' @export
#'
#' @examples
#'
#' #--------------------------------------
#' # Two-sided one-sample test
#' # theta = 0.5
#' # alpha = 0.05, beta = 0.2
#'
#' n <- size.mean(theta = 0.5, sample = "one.sample",
#'                alternative = "two.sided", alpha = 0.05, beta = 0.2)
#'
#' print(n)
#'
#' #--------------------------------------
#' # Two-sided one-sample test
#' # H0: pi = 0.5, H1: pi != 0.5
#' # alpha = 0.05, beta = 0.2, delta = 0.2
#'
#' n <- size.prop(delta = 0.2, pi = 0.5, sample = "one.sample",
#'                alternative = "two.sided", alpha = 0.05, beta = 0.2)
#'
#' print(n)
#'
#' #--------------------------------------
#  # Two-sided test
#' # H0: rho = 0.3, H1: rho != 0.3
#' # alpha = 0.05, beta = 0.2, delta = 0.2
#'
#' n <- size.cor(delta = 0.2, rho = 0.3, alpha = 0.05, beta = 0.2)
#'
#' print(n)
print.size <- function(x, ...) {

  #-----------------------------------------------------------------------------------
  # Main function

  # arithmeetic mean
  if (x$type == "mean") {

    cat("\nSample size determination for the", ifelse(x$spec$sample == "one.sample", "one-sample", "two-sample"), "t-test\n\n")

    if (x$spec$sample == "one.sample") {

        cat(" optimal sample size: n =", ceiling(x$res$n), "\n\n")

      } else {

        cat(" optimal sample size: n =", ceiling(x$res$n), "(in each group) \n\n")

    }

    ###

    # one-sample
    if (x$spec$sample == "one.sample") {

      if (x$spec$alternative == "two.sided") {

        cat("  H0: mu = mu.0  versus  H1: mu != mu.0\n")

      }

      if (x$spec$alternative == "less") {

        cat("  H0: mu >= mu.0  versus  H1: mu < mu.0\n")

      }

      if (x$spec$alternative == "greater") {

        cat("  H0: mu <= mu.0  versus  H1: mu > mu.0\n")

      }

    # two-sample
    } else {

      if (x$spec$alternative == "two.sided") {

        cat("  H0: mu.1 = mu.2  versus  H1: mu.1 != mu.2\n")

      }

      if (x$spec$alternative == "less") {

        cat("  H0: mu.1 >= mu.2  versus  H1: mu.1 < mu.2\n")

      }

      if (x$spec$alternative == "greater") {

        cat("  H0: mu.1 <= mu.2  versus  H1: mu.1 > mu.2\n")

      }

    }

    ###

    cat("  alpha:", x$spec$alpha, " beta:", x$spec$beta, " theta:", x$spec$theta, "\n\n")

  }

  #----------------------------------------

  # proportion
  if (x$type == "prop") {

    if (x$spec$correct == TRUE) {

      cat("\nSample size determination for the", ifelse(x$spec$sample == "one.sample", "one-sample", "two-sample"), "proportion test with continuity correction\n\n")

    } else {

      cat("\nSample size determination for the", ifelse(x$spec$sample == "one.sample", "one-sample", "two-sample"), "proportion test without continuity correction\n\n")

    }

    ###

    if (x$spec$sample == "one.sample") {

      cat(" optimal sample size: n =", ceiling(x$res$n), "\n\n")

    } else {

      cat(" optimal sample size: n =", ceiling(x$res$n), "(in each group) \n\n")

    }

    ###

    # one-sample
    if (x$spec$sample == "one.sample") {

      if (x$spec$alternative == "two.sided") {

        cat("  H0: pi =", x$spec$pi, " versus  H1: pi !=",  x$spec$pi, "\n")

      }

      if (x$spec$alternative == "less") {

        cat("  H0: pi >=", x$spec$pi, " versus  H1: pi <",  x$spec$pi, "\n")

      }

      if (x$spec$alternative == "greater") {

        cat("  H0: pi <=", x$spec$pi, " versus  H1: pi >",  x$spec$pi, "\n")

      }

    # two-sample
    } else {

      if (x$spec$alternative == "two.sided") {

        cat("  H0: pi.1 = pi.2  versus  H1: pi.1 != pi.2\n")

      }

      if (x$spec$alternative == "less") {

        cat("  H0: pi.1 >= pi.2  versus  H1: pi.1 < pi.2\n")

      }

      if (x$spec$alternative == "greater") {

        cat("  H0: pi.1 <= pi.2  versus  H1: pi.1 > pi.2\n")

      }

    }

    ###

    cat("  alpha:", x$spec$alpha, " beta:", x$spec$beta, " delta:", x$spec$delta, "\n\n")

  }

  #----------------------------------------

  # correlation coefficient
  if (x$type == "cor") {

    cat("\nSample size determination for Pearson's correlation coefficient\n\n",

        " optimal sample size: n =", ceiling(x$res$n), "\n\n")

    ###

    if (x$spec$alternative == "two.sided") {

      cat("  H0: rho =", x$spec$rho, " versus  H1: rho !=",  x$spec$rho, "\n")

    }

    if (x$spec$alternative == "less") {

      cat("  H0: rho >=", x$spec$rho, " versus  H1: rho <",  x$spec$rho, "\n")

    }

    if (x$spec$alternative == "greater") {

      cat("  H0: rho <=", x$spec$rho, " versus  H1: rho >",  x$spec$rho, "\n")

    }

    ###

    cat("  alpha:", x$spec$alpha, " beta:", x$spec$beta, " delta:", x$spec$delta, "\n\n")

  }

  #-----------------------------------------------------------------------------------

}
