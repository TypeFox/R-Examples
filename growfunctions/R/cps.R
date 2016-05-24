#' Monthly employment counts from 1990 - 2013 from the Current Population Survey
#' 
#' Monthly employment counts published by the U.S. Bureau of Labor Statistics in the 
#' Current Population Survey (CPS) for each of \code{N = 51} states (including the District
#' of Columbia).  This dataset covers \code{T = 278} months from \emph{1990} the first two
#' months of \emph{2013}.  The data include a {N x T} matrix, \code{y_raw}, of raw employment counts, 
#' as well as set of standardized values, \code{y}, where the standardization is done within state.  
#' The standardized data matrix is used in our \code{\link{gpdpgrow}} and \code{\link{gmrfdpgrow}}
#' estimating functions because the standardization facilitates comparisons of the time-series across
#' states.
#' 
#' \itemize{
#'   \item y. An \code{(N = 51) x (T = 278)} matrix of standardized employment count 
#'         estimates for \code{N = 51} states for \code{T = 278} months, beginning in 1990. The counts
#'         are standardized to (0,1) for each state series
#'   \item y_raw. An \code{N x T} matrix of estimated monthly employment counts for \code{N = 51}
#'         states.
#'   \item st. Two-digit labels for each of the \code{N} states in the order presented in 
#'         \code{y} and \code{y_raw}.
#'   \item dte. A \code{Date} vector of length {T} that presents the set of dates (in \code{y-m-d} 
#'         format) associated to the \code{T} time points presented in \code{y} and \code{y_raw}.
#'   \item yr. A number vector listing sequence of years, \emph{1990 - 2013} included in the data set.
#'   \item yr_label. A numerical vector of length \code{T = 278} with year labels for each 
#'         monthly employment count in the \code{cps} data set. 
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name cps
#' @usage cps
#' @format A list object of 5 objects supporting a data matrix of N = 51 state time series for
#'        T = 278 months.
NULL

 