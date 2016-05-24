#' Sample Size Calculations for Two-Sample Microarray Experiments with 
#' Differing Mean Expressions but fixed Standard Deviations Among Genes
#' 
#' For given desired power, controlled false discovery rate, 
#' and user-specified proportions of non-differentially expressed genes, 
#' \code{ssize.twoSampVaryDelta} calculates appropriate sample sizes for 
#' two-sample microarray experiments in which the differences between mean 
#' treatment expression levels (\emph{delta.g} for gene \emph{g}) 
#' vary among genes. 
#' A plot of power versus sample size is generated.
#' 
#' @importFrom stats pt optimize 
#' @importFrom graphics plot lines points abline title mtext legend
#' 
#' @param deltaMean location (mean) parameter of normal distribution 
#'                  followed by each \emph{delta.g}.
#' @param deltaSE scale (standard deviation) parameter of normal distribution 
#'                followed by each \emph{delta.g}.
#' @param sigma the common standard deviation of expressions for all genes.
#' @param fdr the false discovery rate to be controlled.
#' @param power the desired power to be achieved.
#' @param pi0 a vector (or scalar) of proportions of non-differentially 
#'            expressed genes.
#' @param maxN the maximum sample size used for power calculations.
#' @param side options are "two-sided", "upper", or "lower".
#' @param cex.title controls size of chart titles.
#' @param cex.legend controls size of chart legend.
#' 
#' @details Each \emph{delta.g} is assumed to follow a Normal distribution 
#' with mean \code{deltaMean} and standard deviation \code{deltaSE}. 
#' The standard deviations of expressions are assumed identical for all genes.
#'
#' If a vector is input for \code{pi0}, sample size calculations 
#' are performed for each proportion.
#' 
#' @return \item{ssize}{sample sizes (for each treatment) at 
#'                      which desired power is first reached.}
#' @return \item{power}{power calculations with corresponding 
#'                      sample sizes.}
#' @return \item{crit.vals}{critical value calculations with 
#'                          corresponding sample sizes.}
#' 
#' @author Ran Bi \email{biran@@iastate.edu}, 
#'         Peng Liu \email{pliu@@iastate.edu}
#' 
#' @references Liu, P. and Hwang, J. T. G. (2007) Quick calculation for 
#' sample size while controlling false discovery rate with application 
#' to microarray analysis. \emph{Bioinformatics} 23(6): 739-746. 
#' 
#' Orr, M. and Liu, P. (2009) Sample size estimation while controlling 
#' false discovery rate for microarray experiments using ssize.fdr package. 
#' \emph{The R Journal}, 1, 1, May 2009, 47-53. 
#' 
#' @seealso \code{\link{ssize.twoSamp}}, \code{\link{ssize.twoSampVary}}, 
#'          \code{\link{ssize.oneSamp}}, \code{\link{ssize.oneSampVary}}, 
#'          \code{\link{ssize.F}}, \code{\link{ssize.Fvary}}
#' 
#' @examples
#' dm <- 1.2; ds <- 0.1  ## the delta.g's follow a Normal(1.2, 0.1) distribution
#' s <- 1                ## common standard deviation
#' fdr <- 0.05           ## false discovery rate to be controlled
#' pwr <- 0.8            ## desired power
#' pi0 <- c(0.5, 0.8, 0.99) ## proportions of non-differentially expressed genes
#' N <- 35               ## maximum sample size for calculations
#' 
#' size <- ssize.twoSampVaryDelta(deltaMean = dm, deltaSE = ds, sigma = s, 
#'                                fdr = fdr, power = pwr, pi0 = pi0, 
#'                                maxN = N, side = "two-sided")
#' size$ssize                ## first sample size(s) to reach desired power
#' size$power                ## calculated power for each sample size
#' size$crit.vals            ## calculated critical value for each sample size
#' 
#' @export
#' 
ssize.twoSampVaryDelta <- function (deltaMean, deltaSE, sigma, fdr = 0.05, 
                                    power = 0.8, pi0 = 0.95, maxN = 35, 
                                    side = "two-sided", cex.title = 1.15, 
                                    cex.legend = 1){
  
  N <- maxN
  a <- fdr
  p <- pi0
  if (side == "two-sided") {
    TSVaryDelta <- function(c) {
      r <- a * (1 - p)/((1 - a) * p)
      dif <- abs((2 * pt(q = -c, df = 2 * n - 2)/
                    (1 - pt(q = c/sqrt(deltaSE^2 / (sigma^2 * 2/n) +1), 
                            df = 2 * n - 2, 
                            ncp = deltaMean/sqrt(deltaSE^2 + sigma^2 * 2/n))
                     + pt(q = -c/sqrt(deltaSE^2 / (sigma^2 * 2/n) +1), 
                          df = 2 * n - 2, 
                          ncp = deltaMean/sqrt(deltaSE^2 + sigma^2 * 2/n))) 
                  - r))
      return(dif)
    }
  }
  if (side == "upper") {
    TSVaryDelta <- function(c) {
      r <- a * (1 - p)/((1 - a) * p)
      dif <- abs((2 * pt(q = -c, df = 2 * n - 2)/
                    (1 - pt(q = c/sqrt(deltaSE^2 / (sigma^2 * 2/n) +1), 
                            df = 2 * n - 2, 
                            ncp = deltaMean/sqrt(deltaSE^2 + sigma^2 * 2/n)))
                  - r))
      return(dif)
    }
  }
  if (side == "lower") {
    TSVaryDelta <- function(c) {
      r <- a * (1 - p)/((1 - a) * p)
      dif <- abs((2 * pt(q = -c, df = 2 * n - 2)/
                    pt(q = -c/sqrt(deltaSE^2 / (sigma^2 * 2/n) +1), 
                       df = 2 * n - 2, 
                       ncp = deltaMean/sqrt(deltaSE^2 + sigma^2 * 2/n)) - r))
      return(dif)
    }
  }
  
  pwr2 <- NULL
  crit <- NULL
  ssize <- matrix(0, nrow = length(pi0), ncol = 3)
  colnames(ssize) <- c("pi0", "ssize", "power")
  up.start <- 50
  for (i in 1:length(pi0)) {
    p <- pi0[i]
    up <- up.start
    for (n in 3:N) {
      ci <- optimize(f = TSVaryDelta, interval = c(0, up))$min
      up <- ci
      if (abs(ci - up.start) >= 1) {
        if (side == "two-sided") {
          pwr.new <- (1 - pt(q = ci/sqrt(deltaSE^2 / (sigma^2 * 2/n) +1), 
                             df = 2 * n - 2, 
                             ncp = deltaMean/sqrt(deltaSE^2 + sigma^2 * 2/n))
                      + pt(q = -ci/sqrt(deltaSE^2 / (sigma^2 * 2/n) +1), 
                           df = 2 * n - 2, 
                           ncp = deltaMean/sqrt(deltaSE^2 + sigma^2 * 2/n)))
        }
      }
      if (abs(ci - up.start) < 1) {
        pwr.new <- 0
        ci <- NA
      }
      crit <- c(crit, ci)
      pwr2 <- c(pwr2, pwr.new)
      if (pwr2[(i - 1) * (N - 2) + n - 2] >= power & ssize[i, 1] == 0) {
        ssize[i, ] <- c(p, n, pwr2[(i - 1) * (N - 2) + n - 2])
      }
    }
  }
  ssize[, 1] <- pi0
  if (sum(ssize == 0) > 0) {
    warning("Desired power not achieved for at least one pi0")
  }
  ssize[ssize == 0] <- NA
  pwrMatrix <- matrix(c(3:N, pwr2), ncol = length(pi0) + 1, 
                      byrow = FALSE)
  for (i in 1:length(pi0)) {
    if (i == 1) {
      plot(3:N, pwrMatrix[, i + 1], col = i, xlim = c(0, N), ylim = c(0, 1), 
           xlab = "", ylab = "", pch = 16)
      lines(3:N, pwrMatrix[, i + 1], col = i, lty = i)
    }
    if (i != 1) {
      points(3:N, pwrMatrix[, i + 1], col = i, pch = 16)
      lines(3:N, pwrMatrix[, i + 1], col = i, lty = i)
    }
  }
  abline(h = power, lty = 2, lwd = 2)
  abline(v = 0:N, h = 0.1 * (0:10), col = "gray", lty = 3)
  title(xlab = "Sample size (n)", ylab = "Power")
  mtext(bquote(paste("Average power vs. sample size with fdr=", .(fdr), ",")), 
        cex = cex.title, padj = -1.85)
  mtext(bquote(paste(Delta[g], "~N(", .(round(deltaMean, 4)), ",", 
                     .(round(deltaSE, 4)), ") and ", sigma[g], " = ", 
                     .(round(sigma, 4)))), cex = cex.title, padj = -0.1)
  legend(x = N, y = 0, xjust = 1, yjust = 0, col = 1:i, pch = c(16, 16, 16), 
         lty = 1:length(pi0), legend = as.character(pi0), 
         bg = "white", title = expression(pi[0]), cex = cex.legend)
  pwrMatrix <- round(pwrMatrix, 7)
  colnames(pwrMatrix) <- c("n", as.character(pi0))
  critMatrix <- matrix(c(3:N, crit), ncol = length(pi0) + 1, byrow = FALSE)
  colnames(critMatrix) <- c("n", as.character(pi0))
  ret <- NULL
  ret$ssize <- ssize
  ret$power <- pwrMatrix
  ret$crit.vals <- critMatrix
  return(ret)
}