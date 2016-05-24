#' @name mes
#' @docType data
#' @title Mammography Experience Study data
#'
#' @format A \code{data.frame} with
#' \eqn{412} observations (rows)
#' and \eqn{7} variables (columns).
#'
#' @details
#' A subset of data from a study to assess
#' factors associated with women's
#' knowledge of and attitude towards mammography.
#' \cr \cr
#' The observed variable values were modified to
#' protect patient confidentiality.
#' \cr \cr
#' Columns are:
#' \describe{
#'  \item{OBS}{Observation/ identification code.}
#'  \item{ME}{Mammography experience (\code{factor}):
#'   \describe{
#'    \item{0}{never}
#'    \item{1}{within_one_year}
#'    \item{2}{over_one_year_ago}}}
#'  \item{SYMPT}{"You do not need a mammogram unless
#'               you have symptoms" (\code{factor}):
#'   \describe{
#'    \item{1}{stongly_agree}
#'    \item{2}{agree}
#'    \item{3}{disagree}
#'    \item{4}{strongly_disagree}}}
#'  \item{PB}{Perveived benefit of mammography. \cr
#'            This is the sum of five scaled responses,
#'            each on a four point scale. \cr
#'            A low value is indicative of a woman with
#'            strong agreement with the benefits of mammography.}
#'  \item{HIST}{Mother or sister with a history
#'              of breast cancer? (\code{factor}):
#'   \describe{
#'    \item{0}{no}
#'    \item{1}{yes}}}
#'  \item{BSE}{Breast self-exam. \cr
#'             "Has anyone taught you how to
#'              examine your own breasts?" (\code{factor}):
#'   \describe{
#'    \item{0}{no}
#'    \item{1}{yes}}}
#'  \item{DETC}{"How likely is it that a mammogram could find
#'               a new case of breast cancer?" (\code{factor}):
#'   \describe{
#'    \item{1}{not_likely}
#'    \item{2}{somewhat_likely}
#'    \item{3}{very_likely}}}
#' }
#'
#' @keywords datasets
#'
#' @source
#' \href{ftp://ftp.wiley.com/public/sci_tech_med/logistic}{
#'       Wiley FTP}
#' @references
#' \bold{H&L 2nd ed.} Page 265. Table 8.1.
#'
#' Zapka JG, Stoddard A, Maul L, Costanza ME 1991.
#' Interval adherence to mammography screening guidelines.
#' \emph{Medical Care}
#' \bold{29}(8):697--707.
#' \cr
#' JSTOR (free):
#' \cr
#' http://www.jstor.org/stable/3766098
#'
#' Costanza ME, Stoddard AM, Gaw VP, Zapka JG 1992.
#' The risk factors of age and family history and their
#' relationship to screening mammography utilization.
#' \emph{Journal of the American Geriatrics Society}
#' \bold{40}(8):774--778.
#' \href{http://dx.doi.org/10.1111/j.1532-5415.1992.tb01848.x}{
#'       Wiley (paywall)}
#'
#' Zapka JG, Hosmer D, Costanza ME, Harris DR, Stoddard A 1992.
#' Changes in mammography use: economic, need and service factors.
#' \emph{American Journal of Public Health}
#' \bold{82}(10):1345--1351.
#' \href{http://dx.doi.org/10.2105/AJPH.82.10.1345}{AJPH (free)}
#'
NULL
