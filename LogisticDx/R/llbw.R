#' @name llbw
#' @docType data
#' @title Longitudinal Low Birth Weight study data
#'
#' @format A \code{data.frame} with
#' \eqn{200} observations (rows)
#' and \eqn{8} variables (columns).
#'
#' @details A hypothetical data set based on
#' the reference below.
#' \cr
#' The woman age 45 was excluded as an outlier.
#' \cr
#' A hypothetical additional number (\eqn{1} to \eqn{3})
#' of births was generated for each woman,
#' yielding an average of \eqn{2.6} births per woman.
#' \cr \cr
#' This is a subset of the original dataset.
#' \cr \cr
#' Columns are:
#' \describe{
#'  \item{ID}{Identification code.}
#'  \item{BIRTH}{Birth number. \eqn{1} to \eqn{4}.}
#'  \item{SMOKE}{Smoking status during pregnancy (\code{factor}):
#'   \describe{
#'    \item{0}{no}
#'    \item{1}{yes}}}
#' \item{RACE}{Race (\code{factor}):
#'   \describe{
#'    \item{1}{white}
#'    \item{2}{black}
#'    \item{3}{other}}}
#'  \item{AGE}{Age of mother.}
#'  \item{LWT}{Weight of mother (lbs) at last menstrual period.}
#'  \item{BWT}{Birth weight (grams).}
#'  \item{LBW}{Low birth weight? (\code{factor}):
#'   \describe{
#'    \item{0}{BWT > 2500g}
#'    \item{1}{BWT <= 2500g}}}
#' }
#'
#' @keywords datasets
#'
#' @source
#' \href{ftp://ftp.wiley.com/public/sci_tech_med/logistic}{
#'       Wiley FTP}
#' @references
#' \bold{H&L 2nd ed.} Sections 1.6.2 and 8.3.
#'
NULL
