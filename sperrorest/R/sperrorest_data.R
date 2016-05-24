# Describe the data set:

#' J. Muenchow's Ecuador landslide data set
#'
#' Data set created by Jannes Muenchow, University of Erlangen-Nuremberg, Germany.
#' These data should be cited as Muenchow et al. (2012) (see reference below). This
#' publication also contains additional information on data collection and the geomorphology
#' of the area. The data set provded here is (a subset of) the one from the 'natural' part of the RBSF area
#' and corresponds to landslide distribution in the year 2000.
#' @name ecuador
#' @docType data
#' @format a \code{data.frame} with point samples of landslide and non-landslide locations
#' in a study area in the Andes of southern Ecuador.
#' @references Muenchow, J., Brenning, A., Richter, M., 2012. Geomorphic process rates of landslides along a 
#' humidity gradient in the tropical Andes. Geomorphology, 139-140: 271-284.
#'
#' Brenning, A., 2005. Spatial prediction models for landslide hazards: review, comparison and evaluation. 
#' Natural Hazards and Earth System Sciences, 5(6): 853-862.
#'
#' @examples 
#' data(ecuador)
#' str(ecuador)
#' library(rpart)
#' ctrl = rpart.control(cp = 0.02)
#' fit = rpart(slides ~ dem + slope + hcurv + vcurv + 
#'    log.carea + cslope, 
#'   data = ecuador, control = ctrl)
#' par(xpd = TRUE)
#' plot(fit, compress = TRUE, main = "Muenchow's landslide data set")
#' text(fit, use.n = TRUE)
NULL
