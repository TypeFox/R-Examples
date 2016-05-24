#' A Package for Computating the Pairwise Likelihood Augmented Cox Estimator for Left-Truncated Data.
#'
#' This package provides both lower-level \code{C++} functions (\code{PLAC_TI()}, \code{PLAC_TV()} and \code{PLAC_TvR()}) and an R wrapper function \code{PLAC()} to calculate the pairwise likelihood augmented Cox estimator for left-truncated survival data as proposed by Wu et al. (2015).
#'
#' @section Wrapper Function \code{PLAC()}:
#' This \code{R} wrapper function calls different \code{C++} function depending on the covariate types \code{data} has.
#'
#' @section C++ Functions:
#' The three \code{C++} functions \code{PLAC_TI()}, \code{PLAC_TV()} and \code{PLAC_TvR()} provide a direct interface to the algorithm in case that users need to supply more flexible time-dependent coavriates other than indicator functions.
#'
#' @references Wu, F., Kim, S., Qin, J., Saran, R. and Li, Y. (2016) "A Pairwise-Likelihood Augmented Estimator for the Cox Model Under Left-Truncation." (Submitted to \emph{Journal of American Statistical Association}.)
#' @docType package
#' @name plac-package
NULL
#> NULL
