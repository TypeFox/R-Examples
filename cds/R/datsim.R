#' Simulate Data for a Single Response Style
#' 
#' Simulate data containing a single repsponse style.
#' 
#' @param nr.indv Integer giving the number of individuals required in the sample.
#' @param m The number of items.
#' @param scales The rating scale used for all items.
#' @param err.coeff The standard error used in simulating the truncated normal distribution.
#' @param resp.style A set of cut points across the interval [0, 1] defining the response style
#' transformation.
#' @param true.mu Optional vector of length \code{m} giving the true preferences for the items.
#' @param a Lower boundary of the truncation interval for the simulated true preferences.
#' @param b Upper boundary for the truncation interval for the simulated true preferences.
#' @param plot.graph Logical indicating whether to visualize the response style in a plot.
#' @param use.copula Logical indicating whether to simulate dependent items using a copula.
#' @param reverse.thresh A proportion giving the proportion of item preferences which should be
#' reversed to induce a negative association.
#' @param \dots Additional arguments passed to \code{\link{plot}}.
#' @author Pieter C. Schoonees
#' @references Schoonees, P.C., Velden, M. van de & Groenen, P.J.F. (2013).
#' Constrained Dual Scaling for Detecting Response Styles in Categorical Data.
#' (EI report series EI 2013-10). Rotterdam: Econometric Institute.
#' @keywords multivariate
#' @export datsim
datsim <- function (nr.indv = 100, m = 5, scales = 1:7, err.coeff = 0.1, 
    resp.style = c(-Inf, 1/7, 2/7, 3/7, 4/7, 5/7, 6/7, Inf), 
    true.mu = NULL, a = 0, b = 1, plot.graph = FALSE, use.copula = FALSE, 
    reverse.thresh = 1, ...) {
    qq <- length(scales)
    if (length(resp.style) != qq + 1) 
        stop("Scale and response style have \n\tdifferent lengths")
    if (is.null(true.mu)) {
        true.mu <- runif(m)}
    temp.mat <- matrix(true.mu, nrow = nr.indv, ncol = m, 
        byrow = TRUE)
      sim.mat <- matrix(NA, nrow = nr.indv, ncol = m)
    if(!use.copula){
      for(i in 1:m) sim.mat[, i] <- trRnorm(n = nr.indv, 
           mu = true.mu[i], sd = err.coeff, a = a, b = b)
      true.tau <- NULL
    }
    if(use.copula){
      tmp <- gen.cop(n = nr.indv, true.mu = true.mu, err.coeff = err.coeff, 
                     reverse.thresh = reverse.thresh)
      sim.mat <- tmp$samp
      true.tau <- tmp$tau
      rm(tmp)
    }
    resp.mat <- matrix(as.numeric(cut(sim.mat, breaks = resp.style, 
        labels = scales)), nrow = nr.indv, ncol = m)
    # Add boundaries and rank:
    out <- addbounds(x = resp.mat, q = qq)
    if (plot.graph) {
        plot(c(0, 1), c(0, 1), type = "n", asp = 1, pch = 20, 
            xlab = "Transformed scores", ylab = "Simulated data points", ...)
        points(c(0, scales)/max(scales), c(0, resp.style[-c(1, 
            qq + 1)], 1), type = "b")
        abline(h = true.mu, lty = 3, col = "grey")
        rug(as.numeric(sim.mat), side = 2)
        rug(c(0, scales)/max(scales))
        barplot(table(resp.mat)/sum(table(resp.mat)), ylim = c(0, 
            1), main = "Frequency of ratings")
        abline(h = diff(c(0, resp.style[-c(1, qq + 1)], 1)), 
            lty = 3)
    }
    out <- list(Tmat = out, mu = true.mu, prers = sim.mat, postrs = resp.mat, 
        resp.style = resp.style, m = m, scales = scales, 
        tau = true.tau)
    class(out) <- c("dsdata","list")
    out
}
