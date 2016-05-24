#' Summarized multiple mark-recapture data for Lake Sturgeon.
#'
#' The number of Lake Sturgeon (\emph{Acipenser fulvescens}) caught in multiple samples from Black Lake, MI in 1997.  The caught fish were examined for previous  marks, marked (if previously unmarked), and then returned to the population.
#'
#' @name SturgeonBL
#' @docType data
#' 
#' @format A data frame with 6 observations on the following 4 variables:
#' \describe{
#'   \item{t}{Sample number} 
#'   \item{caught}{Total number of fish caught in the sample}
#'   \item{recaptures}{Number of previously marked fish in the sample}
#'   \item{retmarks}{Number of marked fish (previously and newly marked) returned to the population}
#' }
#' 
#' @section Topic(s):
#'   \itemize{
#'     \item Population Size
#'     \item Abundance
#'     \item Mark-Recapture
#'     \item Capture-Recapture
#'     \item Schnabel
#'     \item Schumacher-Eschmeyer
#' }
#' 
#' @concept Abundance 'Population Size' 'Mark-Recapture' 'Capture-Recapture' Schnabel
#'
#' @source From Table 1 in Baker, E.A. and D.J. Borgeson.  1999. Lake sturgeon abundance and harvest in Black Lake, Michigan, 1975-1999.  North American Journal of Fisheries Management. 19:1080-1088.
#'
#' @keywords datasets
#' 
#' @examples
#' data(SturgeonBL)
#' str(SturgeonBL)
#' head(SturgeonBL)
#'
NULL