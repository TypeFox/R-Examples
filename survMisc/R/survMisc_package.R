#' @name survMisc_package
#' @title Miscellaneous Functions for Survival Analysis
#' 
#' @description
#' \tabular{ll}{
#'  Package: \tab survMisc\cr
#'  Type: \tab Package\cr
#'   Version: \tab 0.5\cr
#'  Date: \tab 2015-07-15\cr
#' License: \tab GPL (>= 2) \cr
#' LazyLoad: \tab yes
#' }
#' A collection of functions for the analysis of survival data. These
#' extend the methods already available in \code{package:survival}.
#' \cr
#' The intent is to generate a workspace for some of the common tasks
#' arising in survival analysis.
#'  \cr \cr
#' There are references in many of the functions to the textbooks:
#' \tabular{cl}{
#'  \bold{K&M} \tab Klein J, Moeschberger M (2003). 
#'                  \emph{Survival Analysis}, 2nd edition. \cr
#'             \tab New York: Springer. 
#'                  \href{http://dx.doi.org/10.1007/b97377}{Springer (paywall)}. \cr
#'  \bold{T&G} \tab Therneau TM, Grambsch PM (2000). 
#'                  \emph{Modeling Survival Data: Extending the Cox Model}. \cr
#'             \tab New York: Springer. 
#'                  \href{http://dx.doi.org/10.1007/978-1-4757-3294-8}{Springer (paywall)}.
#' }
#' 
#' \subsection{Notes for developers}{
#'  \itemize{
#'   \item This package should be regarded as 'in development' until
#' release 1.0, meaning that there may be changes to certain function
#' names and parameters, although I will try to keep this to a minimum.
#' As such it is recommended that other packages do \emph{not} depend on or import from this
#' one until at least version 1.0.
#'  \item Naming tends to follow the \strong{camelCase} convention; 
#'   variables within functions are typically alphanumeric e.g. \code{a1 <- 1}.
#'  }
#' }
#' For bug reports, feature requests or suggestions for improvement,
#' please try to submit to \href{https://github.com/dardisco/survMisc/issues}{github}.
#' Otherwise email me at the address below.
#' 
#' @aliases survMisc
#' @docType package
#' @author Chris Dardis \email{christopherdardis@@gmail.com}
#' 
#' @keywords package
#' @concept survival
#' 
#' @importFrom graphics plot grid title arrows segments points
#' @importFrom grDevices dev.new graphics.off
#' @importFrom utils head tail combn data packageVersion
#' @importFrom stats formula as.formula terms is.empty.model model.matrix model.response printCoefmat predict runif pchisq pnorm qnorm anova
#' 
#' @import knitr
#' @import survival
#' @import ggplot2
#' 
#' @import KMsurv
#' @import km.ci
#' 
#' @importFrom data.table as.data.table data.table set setkey setattr setcolorder setnames ':=' copy rbindlist
#' @importFrom zoo na.locf.default
#' @importFrom grid unit.pmax
#' @importFrom gridExtra grid.arrange
#' 
NULL
## the following are imported for their datasets:
## KMsurv
## km.ci
