#' Observed-Score Linking and Equating
#' 
#' Contains methods for observed-score linking and equating under the
#' single-group, equivalent-groups, and nonequivalent-groups with anchor
#' test(s) designs. Equating types include identity, mean, linear, general
#' linear, equipercentile, circle-arc, and composites of these. Equating
#' methods include synthetic, nominal weights, Tucker, Levine observed score,
#' Levine true score, Braun/Holland, frequency estimation, and chained
#' equating. Plotting and summary methods, and methods for multivariate
#' presmoothing and bootstrap error estimation are also provided.
#' 
#' \tabular{ll}{ Package: \tab equate\cr Version: \tab 2.0-4\cr Date: \tab
#' 2016-4-26\cr Depends: \tab R (>= 3.0.0)\cr License: \tab GPL-3\cr }
#' 
#' Index: \tabular{ll}{
#' ACTmath \tab ACT Mathematics Test Scores\cr
#' KBneat \tab Test Scores Under a NEAT design\cr
#' PISA \tab Programme for International Student Assessment 2009 USA Data\cr
#' descriptives \tab Descriptive Statistics for Frequency Tables\cr
#' percentiles \tab Percentile Ranks and Cumulative Frequencies for Frequency Tables\cr
#' equate \tab Observed Score Linking and Equating\cr
#' composite \tab Composite Linking and Equating\cr
#' bootstrap \tab Bootstrap Equating Error\cr
#' freqtab \tab Frequency Distribution Tables\cr
#' presmoothing \tab Frequency Distribution Presmoothing\cr
#' plot.freqtab \tab Plotting Frequency Distributions\cr
#' plot.equate \tab Plotting Equating Results\cr}
#' 
#' Further information is available in the following vignettes:
#' \tabular{ll}{
#' \code{equatevignette} \tab equate vignette (source, pdf)\cr }
#' The package
#' vignette provides an introduction to linking and equating and includes
#' descriptions of the supported equating methods and examples. The help page
#' for the main function of the package, \code{\link{equate}}, contains
#' additional examples.
#' 
#' @name equate-package
#' @docType package
#' @author Anthony Albano <tony.d.albano@@gmail.com>
#' @importFrom stats AIC BIC as.formula coef fitted formula glm lm median
#' model.matrix poisson reformulate sd stat.anova terms
#' @importFrom graphics layout legend lines matlines par plot plot.default
#' points segments
#' @importFrom grDevices col2rgb rainbow rgb
#' @keywords package
NULL
