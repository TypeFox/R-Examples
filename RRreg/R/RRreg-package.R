#' Correlation and Regression Analyses for Randomized Response Designs
#' 
#' Univariate and multivariate methods for randomized response (RR) survey designs (e.g., Warner, 1965). Univariate estimates of true proportions can be obtained using \code{\link{RRuni}}. RR variables can be used in multivariate analyses for correlations (\code{\link{RRcor}}), as dependent variable in a logistic regression (\code{\link{RRlog}}), or as predictors in a linear regression (\code{\link{RRlin}}). The function \code{\link{RRgen}} generates single RR data sets, whereas \code{\link{RRsimu}} generates and analyzes RR data repeatedly for simulation and bootstrap purposes. An overview of the available RR designs and examples can be found in the package vignette by \code{vignette('RRreg')}.
#'  
#' @details
#' \tabular{ll}{
#' Package: \tab RRreg\cr
#' Type: \tab Package\cr
#' Version: \tab 0.6.0\cr
#' Date: \tab 2015-12-14\cr
#' Depends: \tab R (>= 3.0.0)\cr
#' Imports: \tab parallel, doParallel, foreach, stats, grDevices, graphics, lme4\cr
#' Suggests: \tab knitr\cr
# Encoding: \tab UTF-8\cr
#' License: \tab GPL-2\cr
# LazyLoad: \tab yes\cr
#' URL: \tab \url{http://psycho3.uni-mannheim.de/Home/Research/Software/RRreg}\cr
# Vignette: \tab \url{https://dl.dropboxusercontent.com/u/21456540/RRreg/index.html}\cr
#' }
#' 
#' @aliases RRreg-package RRreg
#' @name RRreg-package
#' @docType package
# @title The RRreg Package
#' @author Daniel W. Heck \email{dheck@@mail.uni-mannheim.de} and Morten Moshagen \email{moshagen@@uni-kassel.de}
#' @keywords package
#' @references Warner, S. L. (1965). Randomized response: A survey technique for eliminating evasive answer bias. \emph{Journal of the American Statistical Association, 60}, 63–69.
NULL

#' Minaret Data
#' 
#' Data by Radukic and Musch (2014)
#' 
#' 
#' The following variables are included:
#' \itemize{
#'   \item \code{age} in years
#'   \item \code{leftRight} political left-right orientation on a scale from -5 to 5
#'   \item \code{rrt} response to RR question (SLD with randomization probabilities \code{p=c(2/12,10/12)})
#'   \item \code{condition} group membership in SLD (either randomization probability \code{p[1]} or \code{p[2]})
#'   \item \code{RRdesign} whether the respondent answered to the RR question (RRdesign=1) or to the direct question (RRdesign=-1)
#'   \item \code{leftRight.c} zero-centered political left-right orientation
#'   \item \code{age.c} zero-centered age
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name minarets
#' @usage data(minarets)
# @format A data frame with 1621 rows and 6 variables
NULL
