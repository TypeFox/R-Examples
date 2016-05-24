#' @name bbdm
#' @docType data
#' @title Benign Breast Disease Matched study data
#'
#' @format A \code{data.frame} with
#' \eqn{200} observations (rows)
#' and \eqn{14} variables (columns).
#'
#' @details The relationship between the use of oral
#' contraceptives and fibrocystic breast disease
#' was examined in a hospital-based case-control study
#' undertaken in New Haven, Connecticut, from 1977 to 1979.
#' \cr \cr
#' This is a subset of the original dataset.
#' \cr \cr
#' Columns are:
#' \describe{
#'  \item{STR}{stratum \eqn{1-50}).}
#'  \item{OBS}{observation within stratum (\code{factor}):
#'   \describe{
#'    \item{1}{Case}
#'    \item{2-4}{Control}}}
#'  \item{AGMT}{Age (years) at interview.}
#'  \item{FNDX}{Final diagnosis (\code{factor}):
#'   \describe{
#'    \item{0}{Control}
#'    \item{1}{Case}}}
#'  \item{HIGD}{Highest grade in school. \eqn{5-20}.}
#'  \item{DEG}{Degree (\code{factor}):
#'   \describe{
#'    \item{0}{none}
#'    \item{1}{high_school}
#'    \item{2}{junior_college}
#'    \item{3}{college}
#'    \item{4}{masters}
#'    \item{5}{doctoral}}}
#'  \item{CHK}{Regular medical checkups? (\code{factor}):
#'   \describe{
#'    \item{1}{Yes}
#'    \item{2}{No}}}
#'  \item{AGP1}{Age (years) at first pregnancy.}
#'  \item{AGMN}{Age (years) at menarche.}
#'  \item{NLV}{Non-live 'births'.
#'             Number of stillbirths, miscarraiges etc.
#'             \eqn{0-7}.}
#'  \item{LIV}{Number of live births. \eqn{0-11}.}
#'  \item{WT}{Weight (lbs) at time of interview.}
#'  \item{AGLP}{Age (years) at last menstrual period.}
#'  \item{MST}{Marital status (\code{factor}):
#'   \describe{
#'    \item{1}{married}
#'    \item{2}{divorced}
#'    \item{3}{separated}
#'    \item{4}{widowed}
#'    \item{5}{never_married}}}
#' }
#'
#' @keywords datasets
#'
#' @source
#' \href{ftp://ftp.wiley.com/public/sci_tech_med/logistic}{
#'       Wiley FTP}
#' @references
#' Pastides H, Kelsey JL, LiVolsi VA, Holford TR,
#' Fischer DB, Goldenberg IS 1983.
#' Oral contraceptive use and fibrocystic breast disease
#' with special reference to its histopathology.
#' \emph{Journal of the National Cancer Institute}
#' \bold{71}(1):5--9.
#' \href{http://jnci.oxfordjournals.org/content/71/1/5.2}{
#'       Oxford (paywall)}
#'
#' Pastides H, Kelsey JL, Holford TR, LiVolsi VA 1985.
#' The epidemiology of fibrocystic breast disease
#' with special reference to its histopathology.
#' \emph{American Journal of Epidemiology}
#' \bold{121}(3):440--447.
#' \href{http://aje.oxfordjournals.org/content/121/3/440}{
#'       Oxford (paywall)}
NULL

