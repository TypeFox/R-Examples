#' Horse mackerel (Trachurus trachurus) off the portuguese coast
#'
#' A dataset containing data collected in two bottom-trawl surveys carried out
#' in the Portuguese coast in 1992 and 1993, and data generated and sampled from
#' it.
#' 
#' Proportions of each age class in the catch (\code{pj}) were calculated using
#' classic Age-Length Keys built from the surveys' data. Mean length at age
#' (\code{lmed}) was given values close to the ones observed in the real data,
#' and its standard deviation (\code{stdv}) was calculated as suggested by
#' Schnute and Fournier (1980).
#' 
#' The two catch matrices \code{N1992} and \code{N1993} were generated from
#' this data, and the \code{F1992} and \code{F1993} were calculated from them.
#' 
#' The \code{otoliths} list contains a total of 1000 length-stratified random
#' samples, extracted from the \code{N1992} matrix. Each of these samples
#' simulate process of the sampling of a small subset of fish from the total
#' catch to analyze its otoliths for age determination. The age data then
#' obtained can then be used either alone or in combination with length data to
#' calculate Age-Length Keys for the population.
#' 
#' The data is presented as a list containing the following items:
#' \itemize{
#'   \item lmed. Mean length at age.
#'   \item stdv. Standard deviation of length at age.
#'   \item pj_1992. Proportions of each age-class in the catch (1992 survey).
#'   \item pj_1993. Proportions of each age-class in the catch (1993 survey).
#'   \item N1992. Number of individuals per length and age class (1992 survey).
#'   \item N1993. Number of individuals per length and age class (1993 survey).
#'   \item F1992. Number of individuals per length (1992 survey).
#'   \item F1993. Number of individuals per length (1993 survey).
#'   \item otoliths. A list of 1000 length-stratified samples taken from N1992.
#' }
#' 
#' @references
#' Schnute, J., Fournier, D. (1980). A New Approach to Length-Frequency
#' Analysis: Growth Structure. Canadian Journal of Fisheries and Aquatic
#' Sciences, 37(9), 1337-1351. DOI: 10.1139/f80-172
#' 
#' @docType data
#' @keywords datasets
#' @format A list of 9 elements
#' @name hom
NULL

