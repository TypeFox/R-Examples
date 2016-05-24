#' Anfis' trained example to use for demonstration
#'
#' The example consist in learning of a bidimentional sinc(x,y) function using
#' a regular grid of 121 points in the domain [-10,10]x[-10,10] and  
#' five independent Normalized Gaussian Membership Function (MF) for each input 
#' (x and y). The training process used the Hybrid off-line Jang's strategy for 
#' 10 epochs.
#' 
#' \describe{
#'  \item{Trainning Set}{\itemize{
#'    \item dim(x)= 121x2, the grid points.
#'    \item dim(y)= 121x1, the sinc(x, y) output.}} 
#'  \item{Architecture}{2 ( 5x5 ) - 25 - 75 ( 75x1 ) - 1, i. e., 2 inputs with 
#'    five MFs in each input, 25 rules and 75 consequents for the single output 
#'    (75x1)}
#'  \item{Last training error}{0.01916307}
#' }
#'
#' @include Anfis-trainSet.R
#' @docType data
#' @format An ANFIS trained object for demonstration. 
#' @source see \code{\link{trainSet}}
#' @name anfis3
#' @family ANFIS
#' @author Cristobal Fresno \email{cfresno@@bdmg.com.ar}, Andrea S. Llera 
#'  \email{ALlera@@leloir.org.ar} and Elmer A. Fernandez 
#'  \email{efernandez@@bdmg.com.ar}
NULL