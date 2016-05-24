#' Comparison Distribution Descriptives
#'
#' calculate a comparison distribution and associated descriptive statisics.
#'
#' This function works within \code{compare} to calculate the comparison
#' distribution of strategy 2 against strategy 1. If \code{trim == TRUE}, a
#' truncated distribution consisting of +/- 2 * interquartile range around the
#' median precision ratio.
#'
#' The victory rate for linear regression is the left tail probability with
#' a cut-off of 1. The victory rate for logistic regression instead taken
#' at a cut-off of 0. For more details, refer to the accompanying PDF document.
#'
#' If \code{output == TRUE} the function will plot a histogram of the distribution
#' and a scatterplot of the precision (SSE or -2 log likelihood) of strategy 2 against
#' strategy 1. Points that lie below the red or blue line represent trials in which
#' strategy 2 outperformed strategy 1 (the victory rate).
#'
#' @importFrom stats glm.fit binomial median IQR
#' @import graphics
#'
#' @param a    a vector containing the sum of square errors or
#'             -2 * log likelihood estimates derived using modelling strategy 1.
#' @param b    a vector containing the sum of square errors or
#'             -2 * log likelihood estimates derived using modelling strategy 2.
#' @param model    type of regression model. Either "linear" or "logistic".
#' @param output   logical. If output is TRUE the function will return two
#'                 graphical representations of the comparison distribution.
#' @param lambda1  a vector or matrix containing the estimated shrinkage factors
#'                 derived using strategy 1
#' @param lambda2  a vector or matrix containing the estimated shrinkage factors
#'                 derived using strategy 2
#' @param trim     logical. If trim is TRUE a "trimmed" comparison distribution
#'                 will be returned, along with a victory rate and median precision
#'                 ratio derived using the trimmed distribution. The trimmed
#'                 distribution only contains precision ratios within a range of
#'                 plus or minus two times the interquartile range around the
#'                 median precision ratio.
#' @param strat1   a list containing the strategy name and strategy-specific
#'                 parameter values.
#' @param strat2   a list containing the strategy name and strategy-specific
#'                 parameter values.
#'
#' @note For details of returned values and examples see \code{\link{compare}}.

compdist <- function(a, b, model, output, lambda1, lambda2, trim, strat1, strat2) {

  if (missing(lambda1)) lambda1 <- "N/A"
  if (missing(lambda2)) lambda2 <- "N/A"
  name1 <- strat1[1]
  name2 <- strat2[1]

  if (model == "linear") {
    distr <- sqrt(b / a)
    VR <- length(subset(distr, distr < 1)) / length(distr)
    pr.med <- round(median(distr), 2)
    pr.quart <- round(IQR(distr), 2)
    if(trim == TRUE) {
      distr.trim <- distr[which(distr > (pr.med - 2 * pr.quart) & distr <
                                  (pr.med + 2 * pr.quart))]
      N.rejected <- length(distr) - length(distr.trim)
      VR.trim <- length(subset(distr.trim, distr.trim < 1)) / length(distr.trim)
      pr.med.trim <- round(median(distr.trim), 2)
    }

  } else if (model == "logistic") {

    distr <- b - a
    # This is -2logLikilihood(b) - -2loglikelihood(a)
    # The victory rate is the proportion of PRs < 0
    # or the left tail with cut-off 0
    VR <- length(subset(distr, distr < 0)) / length(distr)
    pr.med <- round(median(distr, na.rm = TRUE), 2)
    pr.quart <- round(IQR(distr, na.rm = TRUE), 2)
    if (trim == TRUE) {
      distr.trim <- distr[which(distr > (pr.med - 2 * pr.quart) & distr <
                                  (pr.med + 2 * pr.quart))]
      N.rejected <- length(distr) - length(distr.trim)
      VR.trim <- round(length(subset(distr.trim, distr.trim < 0)) / length(distr.trim), 2)
      pr.med.trim <- round(median(distr.trim), 2)
    }
  }

  if(output == TRUE) {
    if (trim == TRUE) {
      hist(distr.trim, breaks = 20, main = "Comparison Distribution (Trimmed)",
           xlab = sprintf("Precision Ratio of %s over %s", name2, name1),
           yaxt = 'n', ylab = NULL, col = "grey")
      } else {
        hist(distr, breaks = 20, main = "Comparison Distribution",
             xlab = sprintf("Precision Ratio of %s over %s", name2, name1),
             yaxt = 'n', ylab = NULL, col = "grey")
      }
    if (model == "linear") {
      abline(v = 1, col = "red", lwd = 2)
      plot(a, b, main = sprintf("Performance of %s over %s", name2, name1),
           xlab = sprintf("%s model SSE", name1),
           ylab = sprintf("%s model SSE", name2))
      abline(0, 1, col = "red")
      } else {
        abline(v = 0, col = "blue", lwd = 2)
        plot(a, b, main = sprintf("Performance of %s over %s", name2, name1),
             xlab = sprintf("%s model -2LogLik", name1),
             ylab = sprintf("%s model -2 LogLik", name2))
        abline(0, 1, col = "blue")
      }
  }


  if(trim == TRUE) {
    comp.out <- c(paste(sprintf("Trimmed Victory Rate of %s over %s =", name2, name1),
                        VR.trim), paste("Trimmed Median Precision Ratio =", pr.med.trim),
                  paste("Precision Ratio IQR =", pr.quart))
    comp.out <- cat(comp.out, sep = "\n")

    comparison <- list(VR.trim = VR.trim, MPR.trim = pr.med.trim, VR = VR,
                       MPR = pr.med, PR.IQR = pr.quart, distribution.trim = distr.trim,
                       N.rejected = N.rejected, distribution = distr, strat1 = a,
                       strat2 = b, shrinkage1 = lambda1, shrinkage2 = lambda2)

    } else {
  comp.out <- c(paste(sprintf("Victory Rate of %s over %s =", name2, name1), VR),
                paste("Median Precision Ratio =", pr.med),
                paste("Precision Ratio IQR =", pr.quart))
  comp.out <- cat(comp.out, sep = "\n")

  comparison <- list(VR = VR, MPR = pr.med, PR.IQR = pr.quart, distribution = distr,
                     strat1 = a, strat2 = b, shrinkage1 = lambda1,
                     shrinkage2 = lambda2)
    }

}
