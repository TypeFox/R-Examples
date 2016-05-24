#' @name SoberanisCruzData
#' @aliases SoberanisCruzData
#' @docType data
#' @title Randomized Response Survey on speeding
#' 
#' @description This data set contains observations from a randomized response survey conducted in a population of 1500 families in a Spanish town 
#' to investigate speeding. 
#' The sample is drawn by cluster sampling by district.
#' The randomized response technique used is the SoberanisCruz model (Soberanis Cruz et al., 2008) with parameter \eqn{p=0.7}.
#' The innocuous question is: Is your car medium/high quality? with \eqn{\alpha=0.5}.
#' 
#' @format
#' A data frame containing 290 observations from a population of \eqn{N=1500} families divided into twenty cluster. 
#' The variables are:
#' \itemize{
#'     \item ID: Survey ID
#'     \item CL: Cluster ID
#'     \item z: The randomized response to the question: Do you often break the speed limit?   
#'     \item Pi: first-order inclusion probabilities  
#' }
#' 
#' @usage SoberanisCruzData
#' 
#' @examples data(SoberanisCruzData)
#' 
#' @keywords datasets
#' 
#' @references Soberanis Cruz, V., Ramírez Valverde, G., Pérez Elizalde, S., González Cossio, F. (2008).
#' \emph{Muestreo de respuestas aleatorizadas en poblaciones finitas: Un enfoque unificador.} 
#'  Agrociencia Vol. 42 Núm. 5 537-549.
#'  
#' @seealso \code{\link{SoberanisCruz}}
NULL