#' Simulated data from Schafer and Galagate (2015)
#'
#' Simulated data used in the paper "Causal inference with a continuous
#' treatment and outcome: alternative estimators for parametric dose-response
#' models".
#'
#' A dataset containing sim_data.
#'
#' @format A data frame with 1000 rows and 20 variables:
#'
#' @return
#' \code{(A.1, A.2, A.3, A.4, A.5, A.6, A.7, A.8)} are the true measured covariates.
#'
#'
#' \code{(B.1, B.2, B.3, B.4, B.5, B.6, B.7, B.8)} are the transformed covariates.
#'   \item{T}{treatment}
#'   \item{Theta.1}{unit level intercept}
#'   \item{Theta.2}{unit level slope}
#'   \item{Y}{outcome}
#'
#'
#' @source use the \code{draw_sample} function
#'
##' @references Schafer, J.L., Galagate, D.L. (2015).  Causal inference with a
##' continuous treatment and outcome: alternative estimators for parametric
##' dose-response models. \emph{Manuscript in preparation}.
#'
#' @examples
#'
#' ## Example from Schafer (2015).
#' data(sim_data)
#' head(sim_data)
#'
#'
#' @usage
#' data(sim_data)
#'
#'
#'
"sim_data"
