#' Deep Vein Thrombosis Study
#'
#' A data frame containing time to recurrence of thrombosis and several potential
#' prognostic factors measured at baseline for 929 individuals with deep vein
#' thrombosis or unprovoked pulmonary embolism. 147 events of recurrence were
#' observed during a median follow-up time of 37.8 months.
#'
#' @note The data are a modified and partly simulated version of the data set used
#'     by Eichinger et al. (2010) and are available under a GPL-2 license.
#'
#' @format  The data frame contains observations of 929 individuals and the
#'     following variables:
#'  \describe{
#'      \item{pnr}{ patient number. }
#'      \item{time}{ time to recurrence of thrombosis or end of study in months. }
#'      \item{status}{ = 1 recurrence of thrombosis. }
#'      \item{sex}{ gender. }
#'      \item{fiimut}{ factor II G20210A mutation. }
#'      \item{fvleid}{ factor V Leiden mutation. }
#'      \item{log2ddim}{ log2-transformed D-dimer. }
#'      \item{bmi}{ body mass index. }
#'      \item{durther}{ duration of anticoagulation therapy. }
#'      \item{age}{ age in years. }
#'      \item{loc}{ location of first thrombosis: pulmonary embolism (PE), distal,
#'            or proximal deep vein \cr thrombosis. }
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
#'  data("deepvein")
#'  summary(deepvein)
"deepvein"
