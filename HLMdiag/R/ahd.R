#' Methylprednisolone data
#' 
#' Data from a longitudinal study exmaining the effectiveness of Methylprednisolone
#' as a treatment for patients with severe alcoholic hepatitis. Subjects were
#' randomly assigned to a treatment (31 received a placebo, 35 received the 
#' treatment) and serum bilirubin was measures each week for four weeks.
#' 
#' @usage data(ahd)
#' @format A data frame with 330 observations on the following 5 variables:
#' \describe{
#' \item{treatment}{The treatment a  subject received - a factor. Levels are
#' \code{placebo} and \code{treated}.}
#' \item{subject}{Subject ID - a factor.}
#' \item{week}{Week of the study (0--4) - the time variable.}
#' \item{sbvalue}{Serum bilirubin level (in \eqn{\mu}mol/L). }
#' \item{baseline}{A subject's serum bilirubin level at week 0.}
#' }
#' @name ahd
#' @docType data
#' @keywords datasets
#' @source
#'  Vonesh, E. F. and Chinchilli, V. M. (1997) \emph{Linear and Nonlinear Models for the 
#'  Analysis of Repeated Measurements}. Marcel Dekker, New York.
#' @references
#' Carithers, R. L., Herlong, H. F., Diehl, A. M., Shaw, E. W., Combes, B., 
#' Fallon, H. J. & Maddrey, W. C. (1989) Methylprednisolone therapy in 
#' patients with severe alcoholic hepatitis. \emph{Annals of Internal Medicine}, 
#' \bold{110}(9), 685--690.
NULL