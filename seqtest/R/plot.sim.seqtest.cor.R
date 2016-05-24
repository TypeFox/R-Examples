#' Plot sim.seqtest
#'
#' This function plots the \code{sim.seqtest.cor} object
#'
#' @param x           \code{sim.seqtest.cor} object.
#' @param plot.lines  plot lines connecting points withe the x- and y-axis.
#' @param plot.nom    plot line at the nominal alpha.
#' @param ylim        the y limits of the plot.
#' @param type        what type of plot should be drawn (\code{"p"} for points,
#'                    \code{"l"} for lines and \code{"b"} for both).
#' @param pch         plotting character.
#' @param lty         line type.
#' @param lwd         line width.
#' @param ...         further arguments passed to or from other methods.
#'
#' @author
#' Takuya Yanagida \email{takuya.yanagida@@univie.ac.at}
#'
#' @seealso
#' \code{\link{sim.seqtest.cor}}, \code{\link{seqtest.cor}}
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
#' sim.obj.1 <- sim.seqtest.cor(rho.sim = 0.3, k = seq(4, 16, by = 1), rho = 0.3,
#'                              alternative = "greater",
#'                              delta = 0.25, alpha = 0.05, beta = 0.05,
#'                              runs = 10000)
#'
#' plot(sim.obj.1)
#'
#' # Step 2: Determine the optimal nominal type-II-risk based on
#' #         the optimal size of subsamples (k) from step 1
#'
#' sim.obj.2 <- sim.seqtest.cor(rho.sim = 0.55, k = 16, rho = 0.3,
#'                              alternative = "greater",
#'                              delta = 0.25, alpha = 0.05, beta = seq(0.05, 0.15, by = 0.01),
#'                              runs = 10000)
#'
#' plot(sim.obj.2)
#' }
plot.sim.seqtest.cor <- function(x, plot.lines = TRUE, plot.nom = TRUE, ylim = NULL,
                                 type = "b", pch = 19, lty = 1, lwd = 1, ...) {

  #--------------------------------------------------------------------------------------------------------#
  # Input Check

  if (length(x$spec$k) == 1 &  length(x$spec$beta) == 1) {

    stop("Object x contains result of only one k and one beta")

  }

  #--------------------------------------------------------------------------------------------------------#
  # Main function

  # length(k) > 1
  if (length(x$spec$k) > 1) {

    k <- x$spec$k
    alpha.emp <- x$res$alpha.emp

    if (is.null(ylim)) { ylim <- c(0, 0.2) }

    plot(k, alpha.emp, type = type, lty = lty, pch = pch, lwd = lwd,
         xlab = "", ylab = "", ylim = ylim, axes = FALSE)

    axis(1, at = k)
    axis(2, at = seq(ylim[1], ylim[2], by = 0.025))

    mtext("Number of Observations in each Sub-Sample (k)", side = 1, line = 2.25)
    mtext(expression(paste("Estimated Empirical Type-I-Risk (",  alpha[emp], ")")), side = 2, line = 2.25)

    box()

    ###

    if (plot.lines == TRUE) {

      for (i in 1:length(k)) {

        lines(c(k[i], k[i]), c(0, alpha.emp[i]), lty = 2, col = "gray70")
        lines(c(k[i], 0), c(alpha.emp[i], alpha.emp[i]), lty = 2, col = "gray70")

      }

      points(k, alpha.emp, pch = pch, type = type)

    }

    ###

    if (plot.nom == TRUE) {

      alpha.nom <- x$spec$alpha

      lines(c(0, max(k) + 1), c(alpha.nom, alpha.nom), col = "red2")

    }

  # length(beta) > 1
  } else {

    beta <- x$spec$beta
    beta.emp <- x$res$beta.emp

    if (is.null(ylim)) { ylim <- c(0, 0.30) }

    plot(beta, beta.emp, type = type, lty = lty, pch = pch, lwd = lwd,
         xlab = "", ylab = "", ylim = ylim, axes = FALSE)

    axis(1, at = beta)
    axis(2, at = seq(ylim[1], ylim[2], by = 0.025))

    mtext(expression(paste("Nominal Type-II-Risk (",  beta[nom], ")")), side = 1, line = 2.5)
    mtext(expression(paste("Estimated Empirical Type-II-Risk ( ",  beta[emp],")")), side = 2, line = 2.25)

    box()

    ###

    if (plot.lines == TRUE) {

      for (i in 1:length(beta)) {

        lines(c(beta[i], beta[i]), c(0, beta.emp[i]), lty = 2, col = "gray70")
        lines(c(beta[i], 0), c(beta.emp[i], beta.emp[i]), lty = 2, col = "gray70")

      }

      points(beta, beta.emp, pch = pch, type = type)

    }

  }

}
