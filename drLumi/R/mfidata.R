#' Median Fluorescence Intensity data
#'
#' @name mfidata
#' @docType data
#' 
#' @description  Dataset with the median fluorescence intensity values 
#' for three 96-wells plates 
#' (from the same experiment), with 30 analytes information per plate 
#' (multiplex experiment).  
#' 
#' @format A \code{data.frame} with the median fluorescence intensity values of 
#' three differents plates based on 30 analytes perplate. 
#' The variables are:
#' \itemize{
#' \item{\code{plate}}{: plate/batch information}
#' \item{\code{well}}{: well position in plate based on 96 wells plate}
#' \item{\code{analyte}}{: analyte (30 per well)}
#' \item{\code{sample}}{: type of sample (blank, standard, positive 
#' control or sample to be analyzed)}
#' \item{\code{mfi}}{: median fluorescence intensity}
#' }
#' 
#' 
#' @keywords data
NULL

