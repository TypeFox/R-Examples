#' Plot the posterior distributions of the focal parameters from a VM model
#'
#' This function plots the univariate and bivariate (if applicable) distributions
#' of the focal (alpha) parameters from a Variability Model where the variability
#' is used as a predictor in a second-stage model.  The latent variability estimates are
#' referred to as \dQuote{Sigma} and, if used, the latent intercepts are referred
#' to as \dQuote{U}.
#'
#' @param alpha Results from running \code{varian} and \code{extract}ing the
#'   results.
#' @param useU Logical indicating whether to plot the latent intercepts
#'   (defaults to \code{TRUE}).  Must set to \code{FALSE} if not available.
#' @param plot Logical whether to plot the results or just return the grob
#'   for the plots.  Defaults to \code{TRUE}.
#' @param digits Integer indicating how many digits should be used
#'   for displaying p-values
#' @param \dots Additional arguments (not currently used)
#' @return A list containing the \code{Combined} and the \code{Individual} plot objects.
#' @author Joshua F. Wiley <josh@@elkhartgroup.com>
#' @importFrom grid grid.draw
#' @export
#' @keywords hplot
#' @examples
#' # Using made up data because the real models take a long time to run
#' set.seed(1234) # make reproducible
#' vmp_plot(matrix(rnorm(1000), ncol = 2))
vmp_plot <- function(alpha, useU = TRUE, plot = TRUE, digits = 3, ...) {

  alpha <- as.data.frame(alpha)
  n <- ncol(alpha)

  stopifnot(n == 1 | n == 2)

  colnames(alpha) <- c("Est_Sigma", "Est_U")[1:n]

  p.sigma <- ggplot(alpha, aes_string("Est_Sigma")) +
    geom_histogram(fill = 'white', colour = 'black',
      binwidth = diff(range(alpha$Est_Sigma, na.rm=TRUE))/50,
      position = "identity") +
    theme_classic()

  sig.sigma <- empirical_pvalue(alpha$Est_Sigma)

  sig.dat <- data.frame(X = c(0, 0),
    Count = sig.sigma[1:2],
    Level = names(sig.sigma)[1:2],
    Pvalue = c(paste0("P = ", format.pval(sig.sigma["p-value"], digits = digits)), ""),
    stringsAsFactors = FALSE)

  if (!useU) {
    p.sig <- ggplot(sig.dat, aes_string("X", "Count", fill = "Level")) +
      geom_bar(stat = 'identity', position = 'stack') +
      scale_fill_manual(values = c("<= 0" = 'grey80', "> 0" = 'grey30')) +
      scale_x_continuous("", breaks = 0, labels = c("Est_Sigma")) +
      geom_text(aes_string("X", "1", label = "Pvalue"), vjust = 0) +
      theme_classic()

    graphs <- list(p.sigma, p.sig)

  } else if (useU) {
    p.u <- ggplot(alpha, aes_string("Est_U")) +
      geom_histogram(fill = 'white', colour = 'black',
        binwidth = diff(range(alpha$Est_U, na.rm=TRUE))/50,
        position = "identity") +
      theme_classic()

    p.joint <- ggplot(alpha, aes_string("Est_Sigma", "Est_U")) +
      geom_point(alpha = .25) + theme_classic()

    sig.u <- empirical_pvalue(alpha$Est_U)

    sig.dat <- rbind(sig.dat, data.frame(X = c(1, 1),
      Count = sig.u[1:2],
      Level = names(sig.sigma)[1:2],
      Pvalue = c(paste0("P = ", format.pval(sig.u["p-value"], digits = digits)), ""),
      stringsAsFactors = FALSE))

    p.sig <- ggplot(sig.dat, aes_string("X", "Count", fill = "Level")) +
      geom_bar(stat = 'identity', position = 'stack') +
      scale_fill_manual(values = c("<= 0" = 'grey80', "> 0" = 'grey30')) +
      scale_x_continuous("", breaks = 0:1, labels = c("Est_Sigma", "Est_U")) +
      geom_text(aes_string("X", "1", label = "Pvalue"), vjust = 0) +
      theme_classic()

    graphs <- list(p.sigma, p.u, p.joint, p.sig)
  }

  p.out <- do.call(arrangeGrob, c(graphs, ncol = 2))

  if (plot) grid.draw(p.out)

  invisible(list(Combined = p.out, Individual = graphs))

}

#' Plot diagnostics from a VM model
#'
#' This function plots a variety of diagnostics from a Variability Model.
#' These include a histogram of the Rhat values (so-called percent scale reduction
#' factors).  An Rhat value of 1 indicates that no reduction in the variability of
#' the estimates is possible from running the chain longer.  Values below 1.10 or 1.05
#' are typically considered indicative of convergence, with higher values indicating
#' the model did not converge and should be changed or run longer.
#' A histogram of the effective sample size indicates for every parameter estimated how
#' many effective posterior samples are available for inference.  Low values may indicate
#' high autocorrelation in the samples and may be a sign of failure to converge.
#' The maximum possible will be the total iterations available.
#' Histograms of the posterior medians for the latent variability and intercept estimates
#' are also shown.
#'
#' @param object Results from running \code{varian}.
#' @param plot Logical whether to plot the results or just return the grob
#'   for the plots.  Defaults to \code{TRUE}.
#' @param \dots Additional arguments not currently used
#' @return A graphical object
#' @author Joshua F. Wiley <josh@@elkhartgroup.com>
#' @export
#' @keywords hplot
#' @examples
#' # Make Me!
vm_diagnostics <- function(object, plot=TRUE, ...) {
  if (inherits(object, "vm")) {
    object <- object$results
  }

  res.s <- as.data.frame(summary(object)$summary)

  est <- extract(object, permute=TRUE)

  p.rhat <- ggplot(res.s, aes_string("Rhat")) +
    geom_histogram(fill = 'white', colour = 'black',
                   binwidth = diff(range(res.s$Rhat))/50,
                   position = "identity") +
    labs(x = "Rhat for all parameters") +
    theme_classic()

  p.neff <- ggplot(res.s, aes_string("n_eff")) +
    geom_histogram(fill = 'white', colour = 'black',
                   binwidth = diff(range(res.s$n_eff))/50,
                   position = "identity") +
    labs(x = "N_Effective for all parameters") +
    theme_classic()

  sigma <- as.data.frame(t(apply(est$Sigma_V, 2, quantile, probs = c(.025, .5, .975), na.rm=TRUE)))
  colnames(sigma) <- c("LL", "Median", "UL")
  sigma <- sigma[order(sigma[, "Median"]), ]
  sigma$Index <- 1:nrow(sigma)

  U <- as.data.frame(t(apply(est$U, 2, quantile, probs = c(.025, .5, .975), na.rm=TRUE)  ))
  colnames(U) <- c("LL", "Median", "UL")
  U <- U[order(U[, "Median"]), ]
  U$Index <- 1:nrow(U)

  p.sigma.h <- ggplot(sigma, aes_string("Median")) +
    geom_histogram(fill = 'white', colour = 'black',
                   binwidth = diff(range(sigma$Median))/50,
                   position = "identity") +
    labs(x = "Median Est_Sigma") +
    theme_classic()

  p.u.h <- ggplot(U, aes_string("Median")) +
    geom_histogram(fill = 'white', colour = 'black',
                   binwidth = diff(range(U$Median))/50,
                   position = "identity") +
    labs(x = "Median Est_U") +
    theme_classic()

  p.sigma <- ggplot(sigma, aes_string("Index", "Median", ymin = "LL", ymax = "UL")) +
    geom_pointrange() +
    labs(y = "Median + 95% CI for Sigma") +
    theme_classic()

  p.u <- ggplot(U, aes_string("Index", "Median", ymin = "LL", ymax = "UL")) +
    geom_pointrange() +
    labs(y = "Median + 95% CI for U") +
    theme_classic()

  p.diag <- arrangeGrob(
    arrangeGrob(p.rhat, p.neff, ncol = 2),
    arrangeGrob(p.sigma.h, p.u.h, ncol = 2),
    p.sigma,
    p.u, ncol = 1)

  if (plot) {
    grid.draw(p.diag)
  }

  invisible(p.diag)
}
