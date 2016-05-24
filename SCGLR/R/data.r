#' @title Sample dataset of abundance of genera in tropical moist forest
#' @name genus
#' @docType data
#' @author CoForChange project 
# @keywords datasets
#' @note The use of this dataset for publication must make reference to the CoForChange project. 
#' @description Genus gives the abundance of 27 common tree genera in the tropical moist forest 
#' of the Congo-Basin and 40 geo-referenced environmental variables on one thousand 8 by 8 km plots
#'  (observations). Each plot's data was obtained by aggregating the data measured on a variable 
#'  number of previously sampled 0.5 ha sub-plots. Geo-referenced environmental variables were 
#'  used to describe the physical factors as well as vegetation characteristics. 
#'  14 physical factors were used pertaining the description of topography, geology and rainfall 
#'  of each plot. Vegetation is characterized through 16-days enhanced vegetation index (EVI) data.
#' @references S. Gourlet-Fleury et al. (2009--2014) CoForChange project: \url{http://www.coforchange.eu/}
#' @references C. Garcia et al. (2013--2015) CoForTips project: \url{http://www.cofortips.org/}
#' @format
#' \tabular{ll}{
#'    \code{gen1 to gen27} \tab abundance of the 27 common genera.\cr
#'    \code{altitude} \tab above-sea level in meters.\cr
#'    \code{pluvio_yr} \tab mean annual rainfall.\cr
#'    \code{forest} \tab classified into seven classes.\cr
#'    \code{pluvio_1 to pluvio_12} \tab monthly rainfalls.\cr
#'    \code{geology} \tab 5-level geological substrate.\cr
#'    \code{evi_1 to evi_23} \tab 16-days enhanced vegetation indexes.\cr
#'    \code{lon and lat} \tab position of the plot centers.\cr
#'    \code{surface} \tab sampled area.\cr
#' }
#' 
NULL
