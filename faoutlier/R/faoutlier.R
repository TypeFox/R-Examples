#' Influential case detection methods for factor analysis and SEM
#'
#' Implements robust Mahalanobis methods, generalized Cook's distances,
#' likelihood ratio tests, model implied residuals, and various
#' graphical methods to help detect and summarize influential
#' cases that can affect exploratory and confirmatory factor analyses.
#'
#' @name faoutlier
#' @docType package
#' @title Influential case detection methods for FA and SEM
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @import stats MASS parallel lattice mvtnorm graphics sem
#' @importFrom lavaan logLik
#' @importFrom methods is
#' @importFrom utils flush.console tail
#' @keywords package
NULL

#' Description of holzinger data
#'
#' A sample of 100 simulated cases from the infamous Holzinger dataset
#' using 9 variables.
#'
#'
#' @name holzinger
#' @docType data
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords data
NULL

#' Description of holzinger data with 1 outlier
#'
#' A sample of 100 simulated cases from the infamous Holzinger dataset
#' using 9 variables, but with 1 outlier added to the dataset. The first row was replaced
#' by adding 2 to five of the observed variables (odd-numbered items) and subtracting 2 from
#' the other four observed variables (even-numbered items).
#'
#'
#' @name holzinger.outlier
#' @docType data
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords data
NULL
