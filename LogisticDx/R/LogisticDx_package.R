#' Diagnostic Tests for Models with a Binomial Response
#'
#' \tabular{ll}{
#'  Package: \tab LogisticDx \cr
#'  Type: \tab Package \cr
#'   Version: \tab 0.2 \cr
#'  Date: \tab 2015-07-01 \cr
#' License: \tab GPL (>= 2) \cr
#' LazyLoad: \tab yes
#' }
#' Diagnostic tests and plots for GLMs (generalized linear models)
#' with binomial/ binary outcomes, particularly logistic regression.
#' \cr \cr
#' The most commonly used functions are likely to be
#' \code{\link{dx}} (diagnostics),
#' \code{\link{plot.glm}} (diagnostic plots) and
#' \code{\link{gof}} (goodness-of-fit tests).
#' \cr \cr
#' There have been changes to many of the functions between
#' Version 0.1 and 0.2 of this package.
#' \cr \cr
#' The package should be regarded as 'in development' until
#' release 1.0, meaning that there may be changes to certain function
#' names and parameters, although I will try to keep this to a minimum.
#' \cr \cr
#' There are references in many of the functions to the textbook:
#' \cr
#' Hosmer D, Lemeshow S (2003).
#' \emph{Applied logistic regression}, 2nd edition.
#' New York: John Wiley & Sons, Inc.
#' \href{http://dx.doi.org/10.1002/0471722146}{Wiley (paywall)},
#' which is herein referred to as \bold{H&L 2nd ed.}
#' \cr \cr
#' For bug reports, feature requests or suggestions for improvement,
#' please try to submit to
#' \href{https://github.com/dardisco/LogisticDx/issues}{github}.
#' Otherwise, email me at the address below.
#'
#' @title Diagnostic Tests for Models with a Binomial Response
#' @docType package
#' @name logisticDx2-package
#' @aliases logisticDx2
#' @author Chris Dardis \email{christopherdardis@@gmail.com}
#' @keywords package
#' @concept diagnostics
#'
#' @importFrom data.table data.table
#' @importFrom data.table setnames
#' @importFrom data.table setattr
#' @importFrom data.table set
#' @importFrom data.table ':='
#' @importFrom data.table rbindlist
#' @importFrom aod wald.test
#' @importFrom statmod glm.scoretest
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom pROC roc
#' @importFrom pROC plot.roc
#'
#' @importFrom speedglm speedglm
#' @importFrom rms lrm
#'
NULL

