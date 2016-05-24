#' @title Biological data for Black Drum from Virginia waters of the Atlantic Ocean, 2001.
#' 
#' @description Biological data (lengths, weights, ages (from otoliths), and sex) for Black Drum (\emph{Pogonias cromis}) from Virginia waters of the Atlantic Ocean, 2001.
#' 
#' @name BlackDrum2001
#' 
#' @docType data
#' 
#' @format A data frame with 141 observations on the following 9 variables.
#'  \describe{
#'    \item{year}{Year of capture (all 2001)}
#'    \item{agid}{Unique identification number}
#'    \item{spname}{Species name (all \dQuote{Black Drum})}
#'    \item{month}{Month of capture}
#'    \item{day}{Day of capture}
#'    \item{weight}{Weight (lbs) -- most are missing}
#'    \item{tl}{Total length (mm)}
#'    \item{sex}{Sex (\code{female}, \code{male}, and \code{unknown})}
#'    \item{otoage}{Age (yrs; from otoliths)}
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Growth
#'    \item von Bertalanffy
#'    \item Weight-Length
#'  }
#'  
#' @concept Growth 'von Bertalanffy' 'Weight-Length'
#' 
#' @note Used in the \href{http://derekogle.com/IFAR/}{Introductory Fisheries Analyses with R} book.
#' 
#' @source Obtained directly from the Virginia Marine Resources Commission via Hank Liao.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(BlackDrum2001)
#' str(BlackDrum2001)
#' head(BlackDrum2001)
#' plot(tl~otoage,data=BlackDrum2001)
#' 
NULL
