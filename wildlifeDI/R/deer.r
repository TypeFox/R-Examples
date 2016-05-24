#' GPS tracking data of two male deer
#'
#' GPS telemetry data for two male deer during a one-week period in March 2005. The 
#' two deer form a male bachelor group, making them an interesting case study for studying dynamic 
#' interaction patterns. The data are a subset of the data used as a case study in Long \emph{et al.} (2014).
#' 
#' The deer data are stored as a single \code{'ltraj'} object; two bursts contain the fixes for two 
#' individuals (deer37 and deer 38). GPS fixes were attempted at a regular sampling frequency of 15 minutes. 
#' For more information on these data  how the deer data was collected or for citation please see the papers 
#' Webb \emph{et al.} (2009, 2010). 
#'
#' @references
#'  Long, J.A., Nelson, T.A., Webb, S.L., Gee, K.L. (2014) A critical examination
#'  of indices of dynamic interaction for wildlife telemetry studies. \emph{Journal 
#'  of Animal Ecology}, \bold{83}: 1216-1233.\cr \cr
#'  Webb, S.L., Gee, K.L., Demarais, S., Strickland, B.K., DeYoung, R.W.
#'  (2009) Efficacy of a 15-strand high-tensile electric fence to control
#'  white-tailed deer movements. \emph{Wildlife Biology in Practice}, \bold{5}, 45-57.\cr\cr
#'  Webb, S.L., Gee, K.L., Strickland, B.K., Demarais, S., DeYoung, R.W.
#'  (2010) Measuring fine-scale white-tailed deer movements and environmental
#'  influences using GPS collars. \emph{International Journal of Ecology},
#'  \bold{2010}, 1-12.
#'
#' @docType data
#' @keywords datasets
#' @format An \code{ltraj} object with two bursts, representing the two different individual deer:
#' \itemize{
#'  \item Deer no. 37 containing 551 fixes.
#'  \item Deer no. 38 containing 567 fixes.
#'  }
#' @name deer
#' @examples
#' data(deer)
#' deer37 <- deer[1]
#' deer38 <- deer[2]
#' plot(deer37)
#' plot(deer38)
#-----------------------------------
NULL
