#' @name pcs
#' @docType data
#' @title Prostate Cancer Study data
#'
#' @format A \code{data.frame} with
#' \eqn{380} observations (rows)
#' and \eqn{9} variables (columns).
#'
#' @details
#' A subset of data from a study of patient with prostate cancer.
#' Variables measured at the baseline patient exam were used to
#' try to determine whether the tumor had penetrated the
#' prostate capsule.
#' \cr \cr
#' The observed variable values were modified to
#' protect patient confidentiality.
#' \cr \cr
#' Columns are:
#' \describe{
#'  \item{ID}{Identification code.}
#'  \item{CAPSULE}{Tumor penetration of prostatic capsule?
#'                 (\code{factor}):
#'   \describe{
#'    \item{0}{no}
#'    \item{1}{yes}}}
#'  \item{AGE}{Age (years).}
#'  \item{RACE}{Race (\code{factor}):
#'   \describe{
#'    \item{1}{white}
#'    \item{2}{black}}}
#'  \item{DPROS}{Digital rectal exam (\code{factor}):
#'   \describe{
#'    \item{1}{no nodule}
#'    \item{2}{unilobar nodule (left)}
#'    \item{3}{unilobar nodule (right)}
#'    \item{4}{bilobar nodule}}}
#'  \item{DCAPS}{Capsular involvement on rectal exam?
#'                 (\code{factor}):
#'   \describe{
#'    \item{0}{no}
#'    \item{1}{yes}}}
#'  \item{PSA}{Prostate Specific Antigen Value (mg/ml).}
#'  \item{VOL}{Tumor volume (cm3)}
#'  \item{GLEASON}{Gleason score (total).
#'                 Range \eqn{0} to \eqn{10}.}
#' }
#'
#' @keywords datasets
#'
#' @source
#' \href{ftp://ftp.wiley.com/public/sci_tech_med/logistic}{
#'       Wiley FTP}
#' @references
#' \bold{H&L 2nd ed.} Page 25. Section 1.6.3.
#'
NULL
