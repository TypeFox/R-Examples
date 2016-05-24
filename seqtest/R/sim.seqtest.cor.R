#' Simulation of the sequential triangular test for Pearson's correlation coefficient
#'
#' This function performs a statistical simulation for the sequential triangular test
#' for Pearson's correlation coefficient.
#'
#' In order to determine the optimal k, simulation is conducted under the H0 condition, i.e., \code{rho.sim} = \code{rho}.
#' Simulation is carried out for a sequence of k values to seek for the optimal k where the empirical alpha is as close
#' as possible to the nominal alpha.
#' In order to determine optimal beta (with fixed k), simulation is conudcted under the H1 condition,
#' i.e., \code{rho.sim} = \code{rho} + \code{delta} or \code{rho.sim} = \code{rho} - \code{delta}.
#' Simulation is carried out for a sequencen of beta values to seek for the optimal beta where the empirical beta
#' is as close as possible to the nominal beta.
#'
#' In order to specify a one-sided test, argument \code{alternative} has to be used (i.e., two-sided tests are conducted by default).
#' Specifying argument \code{alternative = "less"} conducts the simulation for the null hypothesis, H0: \eqn{\rho} >= \eqn{\rho}.0
#' with the alternative hypothesis, H1: \eqn{\rho} < \eqn{\rho}.0; specifying argument \code{alternative = "greater"} conducts the simluation
#' for the null hypothesis, H0: \eqn{\rho} <= \eqn{\rho}.0 with the alternative hypothesis, H1: \eqn{\rho} > \eqn{\rho}.0.
#'
#' @param rho.sim        simulated population correlation coefficient, \eqn{\rho}.
#' @param k              an integer or a numerical vector indicating the number of observations in each sub-sample.
#' @param rho            a number indicating the correlation coefficient under the null hypothesis, \eqn{\rho}.0.
#' @param alternative    a character string specifying the alternative hypothesis,
#' @param delta          minimum difference to be detected, \eqn{\delta}.
#' @param alpha          type-I-risk, \eqn{\alpha}.
#' @param beta           an integer or a numerical vector indicating the type-II-risk, \eqn{\beta}.
#' @param runs           numer of simulation runs.
#' @param m.x            population mean of simulated vector x.
#' @param sd.x           population standard deviation of simulated vector x.
#' @param m.y            population mean of simulated vector y.
#' @param sd.y           population standard deviation of simulated vector y.
#' @param digits         integer indicating the number of decimal places to be displayed.
#' @param output         logical: if \code{TRUE}, output is shown.
#' @param plot           logical: if \code{TRUE}, plot is shown.
#'
#' @author
#' Takuya Yanagida \email{takuya.yanagida@@univie.ac.at},
#'
#' @seealso
#' \code{\link{seqtest.cor}}, \code{\link{plot.sim.seqtest.cor}}, \code{\link{print.sim.seqtest.cor}}
#'
#' @references
#' Schneider, B., Rasch, D., Kubinger, K. D., & Yanagida, T. (2015).
#' A Sequential triangular test of a correlation coefficient's null-hypothesis: 0 \eqn{< \rho \le \rho}0.
#' \emph{Statistical Papers, 56}, 689-699.
#'
#' @return
#' Returns an object of class \code{sim.seqtest.cor} with following entries:
#'
#' \tabular{ll}{
#'   \code{call}      \tab function call \cr
#'   \code{spec}      \tab specification of function arguments \cr
#'   \code{simres}    \tab list with results (for each k or beta) for each run \cr
#'   \code{res}       \tab data.frame with results, i.e., k, alpha.nom (nominal alpha),
#'                         alpha.emp (estimated empirical alpha), beta.nom (nominal beta),
#'                         beta.emp (empirica beta), p.H0 (proportion decision = H0),
#'                         p.H1 (proportion decision = H1), AVN (average number of V),
#'                         ASN (average number of sample pairs) \cr
#' }
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
#' sim.seqtest.cor(rho.sim = 0.3, k = seq(4, 16, by = 1), rho = 0.3,
#'                 alternative = "greater",
#'                 delta = 0.25, alpha = 0.05, beta = 0.05,
#'                 runs = 10000, plot = TRUE)
#'
#' # Step 2: Determine the optimal nominal type-II-risk based on
#' #         the optimal size of subsamples (k) from step 1
#'
#' sim.seqtest.cor(rho.sim = 0.55, k = 16, rho = 0.3,
#'                 alternative = "greater",
#'                 delta = 0.25, alpha = 0.05, beta = seq(0.05, 0.15, by = 0.01),
#'                 runs = 10000, plot = TRUE)
#' }
sim.seqtest.cor <- function(rho.sim, k, rho,
                            alternative = c("two.sided", "less", "greater"),
                            delta, alpha = 0.05, beta = 0.1,
                            runs = 1000, m.x = 0, sd.x = 1, m.y = 0, sd.y = 1,
                            digits = 3, output = TRUE, plot = FALSE) {

  #-----------------------------------------------------------------------------------
  # Input Check

  if (rho.sim < -1 || rho.sim > 1) {

    stop("Correlation coefficient rho.sim out of bound, specify a values between -1 and 1")

  }

  ###

  if (any(k < 4)) {

    stop("At least k = 4 data pairs are needed")

  }

  ###

  if (rho < -1 || rho > 1) {

    stop("Correlation coefficient rho out of bound, specify a values between -1 and 1")

  }

  ###

  if ((length(k) > 1 && length(beta) > 1) || length(rho.sim) != 1 || length(rho) != 1 || length(alpha) != 1) {

    stop("Vector specification is only allowed for k (with fixed beta) or beta (with fixed k)")

  }

  ###

  if(length(k) > 1 && rho.sim != rho) {

    stop("Determine modified k under the H0 condition, i.e., rho.sim equal rho")

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

  if (any(beta <= 0) || any(beta >= 1)) {

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

      if ((rho - delta) <= 0) {

        stop("Value (rho - delta) out of bound")

      }

    } else {

      if ((rho + delta) >= 1) {

        stop("Value (rho + delta) out of bound")

      }

    }

  }

  ###

  if (length(beta) > 1 && rho.sim == rho && alternative == "two.sided") {

    stop("Determine modified beta under the H1 condition, i.e., rho.sim equal (rho + delta) or (rho - delta)")

  }

  if (length(beta) > 1 && rho.sim == rho && alternative == "greater") {

    stop("Determine modified beta under the H1 condition, i.e., rho.sim equal (rho + delta)")

  }

  if (length(beta) > 1 && rho.sim == rho && alternative == "less") {

    stop("Determine modified beta under the H1 condition, i.e., rho.sim equal (rho - delta)")

  }

  #-----------------------------------------------------------------------------------
  # Header

  l.call <- match.call()

  cat("-----------------------------------------------------------------------------\n")
  cat(" Call:    "); print(l.call)
  cat(" Time:   ", time <- paste(Sys.time()), "\n")
  cat(" R:      ", R.version$version.string, "\n")
  cat(" Package:", pkg.version<- paste0("seqtest version ", packageDescription("seqtest")$Version,
                                        " (", packageDescription("seqtest")$Date, ")"), "\n")
  cat("-----------------------------------------------------------------------------\n")

  #-----------------------------------------------------------------------------------
  # Main function

  if (length(k) == 1 && length(beta) == 1) {

    simres <- internal.sim.seqtest.cor(rho.sim = rho.sim, k = k, rho = rho,
                                       alternative = alternative,
                                       delta = delta, alpha = alpha, beta = beta, runs = runs,
                                       m.x = m.x, sd.x = sd.x, m.y = m.y, sd.y = sd.y)

  } else {

    if (length(k) > 1) {

      j <- 1
      simres <- NULL
      for (i in k) {

        if (i == k[1]) {

          simres[[j]] <- internal.sim.seqtest.cor(rho.sim = rho.sim, k = i, rho = rho,
                                                  alternative = alternative,
                                                  delta = delta, alpha = alpha, beta = beta, runs = runs,
                                                  m.x = m.x, sd.x = sd.x, m.y = m.y, sd.y = sd.y)

          j <- j + 1

        } else {

          simres[[j]] <- internal.sim.seqtest.cor(rho.sim = rho.sim, k = i, rho = rho,
                                                  alternative = alternative,
                                                  delta = delta, alpha = alpha, beta = beta, runs = runs,
                                                  m.x = m.x, sd.x = sd.x, m.y = m.y, sd.y = sd.y)
          j <- j + 1

        }

      }

      names(simres) <- paste("k =", k)

    } else {

      j <- 1
      simres <- NULL
      for (i in beta) {

        if (i == k[1]) {

          simres[[j]] <- internal.sim.seqtest.cor(rho.sim = rho.sim, k = k, rho = rho,
                                                  alternative = alternative,
                                                  delta = delta, alpha = alpha, beta = i, runs = runs,
                                                  m.x = m.x, sd.x = sd.x, m.y = m.y, sd.y = sd.y)
          j <- j + 1

        } else {

          simres[[j]] <- internal.sim.seqtest.cor(rho.sim = rho.sim, k = k, rho = rho,
                                                  alternative = alternative,
                                                  delta = delta, alpha = alpha, beta = i, runs = runs,
                                                  m.x = m.x, sd.x = sd.x, m.y = m.y, sd.y = sd.y)
          j <- j + 1

        }

      }

      names(simres) <- paste("beta =", beta)

    }

  }

  #-----------------------------------------------------------------------------------
  # Return object

  # length(k) == 1 & length(beta) == 1
  if (length(k) == 1 && length(beta) == 1) {

    # H0 condition
    if(rho.sim == rho) {

      empiric <- data.frame(k = k,
                            alpha.nom = alpha,
                            alpha.emp = sum(simres[, "H1"]) / runs,
                            p.H0 = sum(simres[, "H0"]) / runs,
                            p.H1 = sum(simres[, "H1"]) / runs,
                            AVN = mean(simres[, "V.m"]),
                            ASN = mean(simres[, "n.fin"]))

      # H1 condition
    } else {

      empiric <- data.frame(k = k,
                            beta.nom = beta,
                            beta.emp = 1 - sum(simres[, "H1"]) / runs,
                            p.H0 = sum(simres[, "H0"]) / runs,
                            p.H1 = sum(simres[, "H1"]) / runs,
                            AVN = mean(simres[, "V.m"]),
                            ASN = mean(simres[, "n.fin"]))

    }

    # length(k) > 1 | length(beta) > 1
  } else {

    # H0 condition
    if(rho.sim == rho) {

      empiric <- data.frame(cbind(k = k,
                                  alpha.nom = alpha,
                                  alpha.emp = sapply(simres, function(x) sum(x[, "H1"]) / runs),
                                  p.H0 = sapply(simres, function(x) sum(x[, "H0"]) / runs),
                                  p.H1 = sapply(simres, function(x) sum(x[, "H1"]) / runs),
                                  AVN = sapply(simres, function(x) mean(x[, "V.m"])),
                                  ASN = sapply(simres, function(x) mean(x[, "n.fin"]))))
      row.names(empiric) <- NULL

      # H1 condition
    } else {

      empiric <- data.frame(cbind(k = k,
                                  beta.nom = beta,
                                  beta.emp =  1 - sapply(simres, function(x) sum(x[, "H1"]) / runs),
                                  p.H0 = sapply(simres, function(x) sum(x[, "H0"]) / runs),
                                  p.H1 = sapply(simres, function(x) sum(x[, "H1"]) / runs),
                                  AVN = sapply(simres, function(x) mean(x[, "V.m"])),
                                  ASN = sapply(simres, function(x) mean(x[, "n.fin"]))))
      row.names(empiric) <- NULL

    }

  }

  object <- list(call = match.call(),
                 spec = list(rho.sim = rho.sim, k = k, rho = rho,
                             alternative = alternative,
                             delta = delta, alpha = alpha, beta = beta,
                             runs = runs, m.x = m.x, sd.x = sd.x, m.y = m.y, sd.y = sd.y,
                             digits = digits),
                 simres = simres,
                 res = empiric)

  class(object) <- "sim.seqtest.cor"

  #-----------------------------------------------------------------------------------
  # Output

  if (output == TRUE) { print(object) }

  if (plot == TRUE && (length(k) > 1 | length(beta) > 1)) { plot(object) }

  return(invisible(object))

}
