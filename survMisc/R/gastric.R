#' @name gastric
#' @docType data
#' @title gastric cancer trial data
#' @format A \code{data.frame} with \eqn{90} rows (observations) and \eqn{3} columns (variables).
#' @details Data from a trial of locally unresectable gastic cancer.
#' \cr
#' Patients (\eqn{n=45} in each group) were randomized to one of two groups:
#' chemotheapy vs. chemotherapy + radiotherapy.
#' \cr
#' Columns are:
#' \describe{
#'  \item{time}{Time, in days}
#'  \item{event}{Death}
#'  \item{group}{Treatment 
#'   \describe{
#'    \item{0}{chemotherapy}
#'    \item{1}{chemotherapy + radiotherapy}
#'   }
#'  }
#' }
#' 
#' @seealso Examples in \code{\link{comp}}
#' 
#' @source Klein J, Moeschberger. Survival Analysis, 2nd edition. Springer 2003.
#' Example 7.9, pg 224.
#' @references Gastrointestinal Tumor Study Group, 1982.
#' A comparison of combination chemotherapy and
#' combined modality therapy for locally advanced gastric carcinoma.
#' \emph{Cancer}. \bold{49}(9):1771-7. 
#' \href{http://dx.doi.org/10.1002/1097-0142(19820501)49:9<1771::AID-CNCR2820490907>3.0.CO;2-M}{
#'  Wiley (free)}.
#' @references Stablein DM, Koutrouvelis IA, 1985.
#' A two-sample test sensitive to crossing hazards in uncensored and singly censored data.
#' \emph{Biometrics}. \bold{41}(3):643-52.
#' \href{http://www.jstor.org/stable/2531284}{JSTOR}.
#' 
#' @examples
#' data("gastric", package="survMisc", verbose=TRUE)
#' head(gastric)
#' 
NULL
