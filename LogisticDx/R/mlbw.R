#' @name mlbw
#' @docType data
#' @title Matched Low Birth Weight data
#'
#' @format A \code{data.frame} with
#' \eqn{112} observations (rows)
#' and \eqn{9} variables (columns).
#'
#' @details This data was collected as
#' part of a larger study at Bayside Medical Center,
#' Springfield, Massachusetts.
#' It contains information on 56 cases
#' (of low birth weight deliveries) and
#' an equal number of age-matched controls.
#' \cr \cr
#' The observed variable values were modified to
#' protect patient confidentiality.
#' \cr \cr
#' A one-to-one matched set was created from the
#' low birth weight data.
#' For each woman who gave birth to a low birth
#' weight baby, a mother of the same age
#' was randomly selected who
#' did not give birth to a low birth weight baby.
#' For three mothers aged \eqn{< 17}, it was not possible
#' to identify a match.
#' \cr \cr
#' Columns are:
#' \describe{
#'  \item{ID}{Identification code.}
#'  \item{LOW}{Low birth weight? (\code{factor}):
#'   \describe{
#'    \item{0}{BWT > 2500g}
#'    \item{1}{BWT <= 2500g}}}
#'  \item{AGE}{Age of mother.}
#'  \item{LWT}{Weight of mother (lbs) at
#'             last menstrual period.}
#' \item{RACE}{Race (\code{factor}):
#'   \describe{
#'    \item{1}{white}
#'    \item{2}{black}
#'    \item{3}{other}}}
#'  \item{SMOKE}{Smoking status during pregnancy
#'               (\code{factor}):
#'   \describe{
#'    \item{0}{no}
#'    \item{1}{yes}}}
#' \item{PTD}{Pre-term delivery previously?
#'            (\code{factor}):
#'   \describe{
#'    \item{0}{no}
#'    \item{1}{yes}}}
#'  \item{HT}{History of hypertension
#'               (\code{factor}):
#'   \describe{
#'    \item{0}{no}
#'    \item{1}{yes}}}
#'  \item{UI}{History of uterine irritability
#'               (\code{factor}):
#'   \describe{
#'    \item{0}{no}
#'    \item{1}{yes}}}
#' }
#' @seealso \code{\link{lbw}}
#' @keywords datasets
#'
#' @source
#' \href{ftp://ftp.wiley.com/public/sci_tech_med/logistic}{
#'       Wiley FTP}
#' @references
#' \bold{H&L 2nd ed.} Page 230. Section 7.3.
#'
NULL
