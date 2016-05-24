#' Draws a sixway plot of MCMC diagnostics.
#'
#' This function produces a variety of diagnostic plots and statistics for MCMC
#' chains.
#'
#' @param chain A numeric vector, \code{\link{mcmc}} object or \code{\link{mcmc.list}} object 
#' (in which case uses its \code{thin} argument, otherwise assumes thinning = 1), storing the
#' MCMC chain for a chosen parameter.
#' @param name The parameter name. If \code{name = NULL}, the column name of \code{chain} will
#' be used, unless that is also \code{NULL} in which case 'x' is used.
#' @param acf.maxlag Maximum lag at which to calculate the auto-correlation
#' function. \code{acf.maxlag = 100} by default. See \code{\link[stats]{acf}}.
#' @param pacf.maxlag Maximum lag at which to calculate the partial
#' auto-correlation function. \code{pacf.maxlag = 10} by default. See \code{\link[stats]{pacf}}.
#' @param ...  Other graphical parameters (see \code{\link[graphics]{par}} for
#' details).
#'
#' @details A variety of plots and statistics are displayed in an R graphic
#' window, including the following:
#' \itemize{
#' \item a trace plot of the plotted trajectory of
#' an MCMC chain for a model parameter;
#' \item a kernel density plot; kernel density estimates are computed
#' using \code{\link[stats]{density}};
#' \item a plotted autocorrelation function (uses \code{\link[stats]{acf}});
#' \item a plotted partial autocorrelation function (uses \code{\link[stats]{pacf}});
#' \item a plot of the estimated
#' Monte Carlo standard error (\code{\link{MCSE}}) of the posterior estimate of the
#' mean against the number of iterations. As MCMC is a
#' simulation-based approach this induces (Monte Carlo) uncertainty due to the
#' random numbers it uses. This uncertainty reduces with more
#' iterations, and is measured by the MCSE, and so this graph details how long
#' the chain needs to be run to achieve a specific MCSE;
#' \item a box contains two contrasting accuracy diagnostics. The
#' Raftery-Lewis diagnostic (\code{\link[coda]{raftery.diag}}) is a diagnostic
#' based on a particular quantile of the distribution. The diagnostic Nhat is
#' used to estimate the length of Markov chain required to estimate a
#' particular quantile (e.g. the 2.5\% and 97.5\% quantiles) to a given
#' accuracy. The Brooks-Draper diagnostic (\code{\link{BD}}) is a diagnostic based on
#' the mean of the distribution. It is used to estimate the length of Markov
#' chain required to produce a mean estimate to k(=2) significant figures with
#' a given accuracy;
#' \item a box of summary
#' statistics including the posterior mean, sd, mode, quantiles and the
#' effective sample size (ESS) of the chain.
#' }
#'
#' @author Zhang, Z., Charlton, C.M.J., Parker, R.M.A., Leckie, G., and Browne,
#' W.J. (2015) Centre for Multilevel Modelling, University of Bristol.
#'
#' @seealso
#' \code{\link{BD}},\code{\link{MCSE}},\code{\link{density}},\code{\link{acf}},\code{\link{pacf}},\code{\link[coda]{raftery.diag}},\code{\link[coda]{effectiveSize}}
#' @examples
#'
#'
#' \dontrun{
#' library(R2MLwiN)
#' # NOTE: if MLwiN not saved in location R2MLwiN defaults to, specify path via:
#' # options(MLwiN_path = 'path/to/MLwiN vX.XX/')
#' # If using R2MLwiN via WINE, the path may look like this:
#' # options(MLwiN_path = '/home/USERNAME/.wine/drive_c/Program Files (x86)/MLwiN vX.XX/')
#'
#' ## Example: tutorial
#' data(tutorial, package = "R2MLwiN")
#'
#' (mymodel <- runMLwiN(normexam ~ 1 + standlrt + (1 + standlrt | school) + (1 | student),
#'                      estoptions = list(EstM = 1, resi.store.levs = 2), data = tutorial))
#'
#' sixway(mymodel@@chains[, "FP_standlrt", drop = FALSE], "beta_1")
#'
#' }
#'
#' @export
sixway <- function(chain, name = NULL, acf.maxlag = 100, pacf.maxlag = 10, ...) {
  args <- list(...)
  
  isMCMC <- TRUE
  if (!coda::is.mcmc(chain) && !coda::is.mcmc.list(chain)) {
    isMCMC <- FALSE
    chain <- coda::mcmc(chain)
  }
  
  if (coda::nvar(chain) > 1) {
    stop("Cannot plot more than one parameter at a time")
  }
  
  if (coda::is.mcmc.list(chain) && coda::nchain(chain) > 1) {
    warning("Sixway does not fully support multiple chains - diagnostics will be on concatenated chains")
  }
  
  if (length(args) > 0 && "mar" %in% names(args)) {
    mar <- args[["mar"]]
  } else {
    mar <- c(4, 4, 2, 1)/2
  }
  
  if (length(args) > 0 && "mgp" %in% names(args)) {
    mgp <- args[["mgp"]]
  } else {
    mgp <- c(1, 0.25, 0)
  }
  
  if (is.null(name)) {
    if (!is.null(coda::varnames(chain))) {
      name <- coda::varnames(chain)[1]
    } else {
      name <- "x"
    }
  }
  
  dev.new()
  mypar <- par(mar = mar, mgp = mgp, ...)
  on.exit(par(mypar))
  split.screen(figs = c(4, 1))
  split.screen(figs = c(1, 2), screen = 1)
  split.screen(figs = c(1, 2), screen = 2)
  split.screen(figs = c(1, 2), screen = 3)
  split.screen(figs = c(1, 1), screen = 4)
  
  screen(5)
  coda::traceplot(chain, xlab = "Iterations", ylab = "parameter", type = "l", tcl = -0.1, cex.axis = 0.8, main = "")
  
  flatchain <- as.matrix(chain)
  
  screen(6)
  dens <- density(flatchain)
  plot(dens, xlab = "parameter value", ylab = "kernel density", main = "", tcl = -0.1, cex.axis = 0.8)
  
  screen(7)
  aa <- acf(flatchain, acf.maxlag, main = "", mgp = c(1, 0.25, 0), tcl = -0.1, cex.axis = 0.8, ylim = c(0, 1))
  rho <- aa$acf[2]
  
  screen(8)
  pacf(flatchain, pacf.maxlag, main = "", mgp = c(1, 0.25, 0), tcl = -0.1, cex.axis = 0.8, ylim = c(0, 1))
  
  screen(9)
  mcse <- MCSE(flatchain, rho, ll = 0.5, ul = 20)
  plot(mcse[, 1], mcse[, 2], type = "l", xlab = "updates", ylab = "MCSE", tcl = -0.1, cex.axis = 0.8)
  
  RL1 <- raftery.diag(flatchain, q = 0.025, r = 0.005, s = 0.95, converge.eps = 0.001)
  RL2 <- raftery.diag(flatchain, q = 0.975, r = 0.005, s = 0.95, converge.eps = 0.001)
  Ndb <- BD(mean(flatchain), var(flatchain), rho, k = 2, alpha = 0.05)
  
  screen(10)
  plot(1, xlim = c(1, 10), ylim = c(1, 5), type = "n", axes = FALSE, xlab = "", ylab = "", frame.plot = TRUE)
  text(5, 4.8, "Accuracy Diagnostics", cex = 1.2)
  if (RL1$resmatrix[1] == "Error") {
    text(5, 4, paste("RL diagnostic only available after ", RL1$resmatrix[2], " updates.", sep = ""), cex = 0.8)
  } else {
    text(5, 4, paste("Raftery-Lewis (quantile) : Nhat =(", RL1$resmatrix[1, "N"], ",", RL2$resmatrix[1, "N"],
                     ")", sep = ""), cex = 0.8)
  }
  text(5, 3, "when q=(0.025,0.975), r=0.005 and s=0.95", cex = 0.8)
  text(5, 2.1, paste("Brooks-Draper (mean) : Nhat =", Ndb), cex = 0.8)
  text(5, 1.2, "when k=2 sigfigs and alpha=0.05", cex = 0.8)
  screen(11)
  plot(1, xlim = c(1, 22), ylim = c(1, 4), type = "n", axes = FALSE, xlab = "", ylab = "", frame.plot = TRUE)
  text(10, 3.8, "Summary Statistics", cex = 1.2)
  quants <- round(quantile(flatchain, c(0.025, 0.05, 0.5, 0.95, 0.975)), 3)
  text(10, 2.9, paste("param name :", name, "posterior mean =", round(mean(flatchain), 3), "SD = ", round(sd(flatchain),
                                                                                                          3), "mode =", round(dens$x[which.max(dens$y)], 3)), cex = 0.8)
  text(10, 2, paste("quantiles : 2.5% =", quants[1], "5% =", quants[2], "50% =", quants[3], "95% =", quants[4],
                    "97.5% =", quants[5]), cex = 0.8)
  if (isMCMC) {
    text(10, 1.2, paste(coda::niter(chain) * coda::thin(chain), "actual iterations storing every ", paste(coda::thin(chain), "th",
                                                                                              sep = ""), " iteration. Effective Sample Size (ESS) =", round(coda::effectiveSize(chain))), cex = 0.8)
  } else {
    text(10, 1.2, paste(length(chain), "actual iterations. Diagnostics assume storing every 1th iteration. Effective Sample Size (ESS) =",
                        round(coda::effectiveSize(chain))), cex = 0.8)
  }
  close.screen(all.screens = TRUE)
}
