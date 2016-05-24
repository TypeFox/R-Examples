#' @name lbw
#' @docType data
#' @title Low Birth Weight study data
#'
#' @format A \code{data.frame} with
#' \eqn{189} observations (rows)
#' and \eqn{11} variables (columns).
#'
#' @details This data was collected as
#' part of a larger study at Bayside Medical Center,
#' Springfield, Massachusetts.
#' It contains information on 189 births
#' to women that were seen in the obsetetrics clinic.
#' \cr \cr
#' The observed variable values were modified to
#' protect patient confidentiality.
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
#' \item{PTL}{Number of previous premature labors.
#'            0 = none.}
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
#'  \item{FTV}{Number of first trimester physician
#'             visits. 0 = none.}
#'  \item{BWT}{Birth weight (grams).}
#' }
#'
#' @keywords datasets
#'
#' @seealso
#' \code{\link{sig}}
#' \code{\link{OR}}
#'
#' @source
#' \href{ftp://ftp.wiley.com/public/sci_tech_med/logistic}{Wiley FTP}
#' @references
#' \bold{H&L 2nd ed.} Page 24. Section 1.6.2.
#'
NULL
