#' @name uis
#' @docType data
#' @title UMARU IMPACT Study data
#'
#' @format A \code{data.frame} with
#' \eqn{575} observations (rows)
#' and \eqn{9} variables (columns).
#'
#' @details
#' A subset of data from the University of Massachusets Aids
#' Research Unit (UMARU) IMPACT study.
#' \cr
#' This came from two concurrent randomized trials of residential
#' treatement for durg abuse, in order to compare planned
#' durations of admission.
#' \cr
#' Site A randomized 444 participants to compare 3 and 6 month
#' stays in a therapeutic community. They were
#' trained to recognize triggers for relapse and taught skills to
#' cope without using drugs.
#' \cr
#' Site B randomized 184 participants to receive either a 6 or 12
#' month stay in a highly structured communal therapeutic
#' community.
#' \cr \cr
#' This is a subset of the original dataset.
#' \cr \cr
#' Columns are:
#' \describe{
#'  \item{ID}{Identification code.}
#'  \item{AGE}{Age (years).}
#'  \item{BECK}{Beck Depression score on admission.}
#'  \item{IVHX}{IV drug use history (\code{factor}):
#'   \describe{
#'    \item{1}{never}
#'    \item{2}{previous}
#'    \item{3}{current}}}
#'  \item{NDRUGTX}{Number of prior drug treatments.
#'                 Range \eqn{5} to \eqn{20}.}
#'  \item{RACE}{Race (\code{factor}):
#'   \describe{
#'    \item{0}{white}
#'    \item{1}{other}}}
#'  \item{TREAT}{Treatment randomization.
#'               'Short' is 3 months in site A, 6 months in site B.
#'               'Long' is 6 months in site A, 12 months in site B.
#'               (\code{factor}):
#'   \describe{
#'    \item{0}{short}
#'    \item{1}{long}}}
#'  \item{SITE}{Assignment treatment site (\code{factor}):
#'   \describe{
#'    \item{0}{A}
#'    \item{1}{B}}}
#'  \item{DFREE}{Remained drug free for 12 months (\code{factor}):
#'   \describe{
#'    \item{0}{no}
#'    \item{1}{yes}}}
#' }
#'
#' @keywords datasets
#'
#' @seealso
#' \code{\link{dx}}
#' \code{\link{plot.glm}}
#'
#' @source
#' \href{ftp://ftp.wiley.com/public/sci_tech_med/logistic}{
#'       Wiley FTP}
#' @references
#' \bold{H&L 2nd ed.} Page 26. Section 1.6.4.
#'
#' McCusker J, Vickers-Lahti M, Stoddard A, Hindin R, Bigelow C,
#' Zorn M, Garfield F, Frost R, Love C, Lewis B 1995.
#' Fischer DB, Goldenberg IS 1983.
#' The effectiveness of alternative planned durations
#' of residential drug abuse treatment.
#' \emph{American Journal of Public Health}
#' \bold{85}(10):1426--1429.
#' \href{http://dx.doi.org/10.2105/AJPH.85.10.1426}{APHA (free)}
#'
#' McCusker J, Bigelow C, Frost R, Garfield F, Hindin R,
#' Vickers-Lahti M, Lewis B 1997. #'
#' The effects of planned duration of residential drug abuse
#' treatment on recovery and HIV risk behavior.
#' \emph{American Journal of Public Health}
#' \bold{87}(10):1637--1644.
#' \href{http://dx.doi.org/10.2105/AJPH.87.10.1637}{APHA (free)}
#'
#' McCusker J, Bigelow C, Vickers-Lahti M, Spotts D,
#' Garfield F, Frost R 1997.
#' Planned duration of residential drug abuse treatment:
#' efficacy versus effectiveness.
#' \emph{Addiction} \bold{92}(11):1467--1478.
#' \href{http://dx.doi.org/10.1111/j.1360-0443.1997.tb02868.x}{
#'       Wiley (paywall)}
#'
NULL
