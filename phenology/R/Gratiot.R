#' Leatherback nest counts from Gratiot et al. (2006) Figure 1
#' @title Leatherback nest counts from Gratiot et al. (2006) Figure 1
#' @author KWATA ONG - 
#' @docType data
#' @name Gratiot
#' @description Leatherback nest counts from Gratiot et al. (2006) Figure 1. These 
#'              data have been collected by the ONG Kwata in French Guiana.\cr
#'              The data have been obtained from the graph of the publication (see reference).
#' @references Gratiot, N., Gratiot, J., de Thoisy, B. & Kelle, L. 2006. 
#'             Estimation of marine turtles nesting season from incomplete 
#'             data ; statistical adjustment of a sinusoidal function. Animal 
#'             Conservation, 9, 95-102.
#' @keywords datasets
#' @usage Gratiot
#' @examples
#' library(phenology)
#' # Read a file with data
#' \dontrun{
#' Gratiot<-read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", header=FALSE)
#' }
#' data(Gratiot)
#' @format data.frame with the morning date in the first column and the nest counts on the second one.
NULL
