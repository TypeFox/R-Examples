#' Print sim.seqtest
#'
#' This function prints the \code{sim.seqtest.cor} object
#'
#' @param x           \code{sim.seqtest.cor} object.
#' @param ...         further arguments passed to or from other methods.
#'
#' @author
#' Takuya Yanagida \email{takuya.yanagida@@univie.ac.at}
#'
#' @seealso
#' \code{\link{sim.seqtest.cor}}, \code{\link{plot.sim.seqtest.cor}}
#'
#' @references
#' Schneider, B., Rasch, D., Kubinger, K. D., & Yanagida, T. (2015).
#' A Sequential triangular test of a correlation coefficient's null-hypothesis: 0 \eqn{< \rho \le \rho}0.
#' \emph{Statistical Papers, 56}, 689-699.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' #---------------------------------------------
#' # Determine optimal k and nominal type-II-risk
#' # H0: rho <= 0.3, H1: rho > 0.3
#' # alpha = 0.01, beta = 0.05, delta = 0.25
#'
#' # Step 1: Determine the optimal size of subsamples (k)
#'
#' sim.obj <- sim.seqtest.cor(rho.sim = 0.3, k = seq(4, 16, by = 1), rho = 0.3,
#'                            alternative = "greater",
#'                            delta = 0.25, alpha = 0.05, beta = 0.05,
#'                            runs = 10000, output = FALSE)
#'
#' print(sim.obj)
#'
#' # Step 2: Determine the optimal nominal type-II-risk based on
#' #         the optimal size of subsamples (k) from step 1
#'
#' sim.obj <- sim.seqtest.cor(rho.sim = 0.55, k = 16, rho = 0.3,
#'                            alternative = "greater",
#'                            delta = 0.25, alpha = 0.05, beta = seq(0.05, 0.15, by = 0.01),
#'                            runs = 10000, output = FALSE)
#'
#' print(sim.obj)
#' }
print.sim.seqtest.cor <- function(x, ...) {

  #-----------------------------------------------------------------------------------
  # Input Check

  if (!inherits(x, "sim.seqtest.cor")) {

    stop("Object is not a sim.seqtest.cor object!")

  }

  #-----------------------------------------------------------------------------------
  # Main function

  cat("\n Statistical Simulation for the Sequential Triangular Test\n\n")

  if (x$spec$alternative == "two.sided") {

    cat("  H0: rho =", x$spec$rho, " versus  H1: rho !=",  x$spec$rho, "\n\n")

  } else {

    if (x$spec$alternative == "less") {

      cat("  H0: rho >=", x$spec$rho, " versus  H1: rho <",  x$spec$rho, "\n\n")

    } else {

      cat("  H0: rho <=", x$spec$rho, " versus  H1: rho >",  x$spec$rho, "\n\n")

    }

  }

  # length(k) == 1 & length(beta) == 1
  if (length(x$spec$k) == 1 & length(x$spec$beta) == 1) {

    cat("   Nominal type-I-risk (alpha):      ", x$spec$alpha, "\n",
        "  Nominal type-II-risk (beta):      ", x$spec$beta, "\n",
        "  Practical relevant effect (delta):", x$spec$delta, "\n",
        "  n in each sub-sample (k):         ", x$spec$k, "\n\n",

        "  Simulated data based on rho:      ", x$spec$rho.sim, "\n",
        "  Simulation runs:                  ", x$spec$runs, "\n")

    if (x$spec$rho.sim == x$spec$rho) {

      cat("\n   Estimated empirical type-I-risk (alpha):", x$res$alpha.emp, "\n",
          "  Average number of steps (AVN):          ", x$res$AVN, "\n",
          "  Average number of sample pairs (ASN):   ", x$res$ASN, "\n\n")

    } else {

      cat("\n   Estimated empirical type-II-risk (beta):", formatC(x$res$beta.emp, digits = x$spec$digits, format = "f"), "\n",
          "  Average number of steps (AVN):          ", formatC(x$res$AVN, digits = x$spec$digits, format = "f"), "\n",
          "  Average number of sample pairs (ASN):   ", formatC(x$res$ASN, digits = x$spec$digits, format = "f"), "\n\n")
    }

  # length(k) != 1 | length(beta) != 1
  } else {

    # length(k) > 1
    if (length(x$spec$k) > 1) {

      cat("   Nominal type-I-risk (alpha):      ", x$spec$alpha, "\n",
          "  Nominal type-II-risk (beta):      ", x$spec$beta, "\n",
          "  Practical relevant effect (delta):", x$spec$delta, "\n\n",

          "  Simulated data based on rho:      ", x$spec$rho.sim, "\n",
          "  Simulation runs:                  ", x$spec$runs, "\n\n")

      cat("  Estimated empirical type-I-risk (alpha):\n")

      ###

      for (i in x$spec$k) {

          cat(paste0("   k = ", i, ": ", formatC(x$res$p.H1[x$res$k == i], digits = x$spec$digits, format = "f")), "\n")

      }

      cat("\n  Average number of steps (AVN):\n")

      for (i in x$spec$k) {

        cat(paste0("   k = ", i, ": ", formatC(x$res$AVN[x$res$k == i], digits = x$spec$digits, format = "f")), "\n")

      }

      cat("\n  Average number of sample pairs (ASN):\n")

      for (i in x$spec$k) {

        cat(paste0("   k = ", i, ": ", formatC(x$res$ASN[x$res$k == i], digits = x$spec$digits, format = "f")), "\n")

      }

    # length(beta) > 1
    } else {

      cat("   Nominal type-I-risk (alpha):      ", x$spec$alpha, "\n",
          "  Practical relevant effect (delta):", x$spec$delta, "\n",
          "  n in each sub-sample (k):         ", x$spec$k, "\n\n",

          "  Simulated data based on rho:      ", x$spec$rho.sim, "\n",
          "  Simulation runs:                  ", x$spec$runs, "\n\n")

      cat("  Estimated empirical type-II-risk (beta):\n")

      ###

      digits <- max(nchar(x$spec$beta)) - 2

      for (i in x$spec$beta) {

        cat(paste0("   Nominal beta = ", formatC(i, format = "f", digits = digits), ": ",
                   formatC(x$res$beta.emp[x$res$beta.nom == i], digits = x$spec$digits, format = "f")), "\n")

      }

      cat("\n  Average number of steps (AVN):\n")

      for (i in x$spec$beta) {

        cat(paste0("   Nominal beta = ", formatC(i, format = "f", digits = digits), ": ",
                   formatC(x$res$AVN[x$res$beta.nom == i], digits = x$spec$digits, format = "f")), "\n")

      }

      cat("\n  Average number of sample pairs (ASN):\n")

      for (i in x$spec$beta) {

        cat(paste0("   Nominal beta = ", formatC(i, format = "f", digits = digits), ": ",
                   formatC(x$res$ASN[x$res$beta.nom == i], digits = x$spec$digits, format = "f")), "\n")

      }

    }

    cat("\n")

  }

}
