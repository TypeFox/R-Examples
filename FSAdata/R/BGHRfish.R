#' @title Fish information from samples collected from Big Hill Reservoir, KS, 2014.
#' 
#' @description Fish information from samples collected from Big Hill Reservoir, KS, in May, 2014.
#' 
#' @name BGHRfish
#' 
#' @docType data
#' 
#' @format A data frame with 266 observations on the following 6 variables.
#'  \describe{
#'    \item{UID}{Unique sample identification number (see \code{\link{BGHRsample}})}
#'    \item{fishID}{Unique fish identification number}
#'    \item{specCode}{Numeric code for each species (\code{116}=\dQuote{Smallmouth Bass}, \code{118}=\dQuote{Largemouth Bass}, and \code{122}=\dQuote{Bluegill})}
#'    \item{length}{Total length (mm)}
#'    \item{weight}{Weight (g)}
#'    \item{count}{Number of fish sampled of that species and length}
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Data Manipulation
#'  }
#'  
#' @concept 'Data Manipulation'
#' 
#' @note Used in the \href{http://derekogle.com/IFAR/}{Introductory Fisheries Analyses with R} book.
#' 
#' @source Obtained directly from Ben Neely.
#' 
#' @seealso See \code{\link{BGHRfish}} for individual fish collected in these samples.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(BGHRfish)
#' str(BGHRfish)
#' head(BGHRfish)
#' 
NULL
