#' German Breast Cancer Study Group
#'
#' A data frame containing the observations from the GBSG study.
#'
#' @format This data frame contains the observations of 686 women:
#' \describe{
#'     \item{id}{patient id.}
#'     \item{htreat}{hormonal therapy, a factor at two levels \code{0} (no) and \code{1} (yes).}
#'     \item{age}{of the patients in years.}
#'     \item{menostat}{menopausal status, a factor at two levels \code{1} (premenopausal) and \code{2} (postmenopausal).}
#'     \item{tumsize}{tumor size (in mm).}
#'     \item{tumgrad}{tumor grade, a ordered factor at levels \code{1 < 2 < 3}.}
#'     \item{posnodal}{number of positive nodes.}
#'     \item{prm}{progesterone receptor (in fmol).}
#'     \item{esm}{estrogen receptor (in fmol).}
#'     \item{rfst}{recurrence free survival time (in days).}
#'     \item{cens}{censoring indicator (\code{0} censored, \code{1} event).}
#' }
#'
#' @references M. Schumacher, G. Basert, H. Bojar,  K. Huebner, M. Olschewski,
#'     W. Sauerbrei, C. Schmoor, C. Beyerle, R.L.A. Neumann and H.F. Rauschecker
#'     for the German Breast Cancer Study Group (1994).
#'     Randomized \eqn{2 \times 2} trial evaluating hormonal treatment
#'     and the duration of chemotherapy in node-positive breast cancer patients.
#'     \emph{Journal of Clinical Oncology}, \bold{12}, 2086--2093.\cr
#'     W. Sauerbrei and P. Royston (1999). Building multivariable prognostic
#'     and diagnostic models: transformation of the predictors by using
#'     fractional polynomials. \emph{Journal of the Royal Statistics Society
#'       Series A}, Volume \bold{162}(1), 71--94.
#'
#' @keywords datasets
#'
#' @examples
#'  data("GBSG")
#'  summary(GBSG)
"GBSG"
