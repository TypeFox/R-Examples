#' @title Summarized multiple mark-recapture data for Redear Sunfish.
#'
#' @description The number of Redear Sunfish (\emph{Lepomis microlophus}) caught in multiple samples from Gordy Lake, IN.  The caught fish were examined for previous marks, marked (if previously unmarked), and then returned to the population.
#'
#' @name SunfishIN
#' 
#' @docType data
#' 
#' @format A data frame with 6 observations on the following 4 variables:
#'  \describe{
#'    \item{t}{Sample number} 
#'    \item{caught}{Total number of fish caught in the sample}
#'    \item{recaps}{Number of previously marked fish in the sample}
#'    \item{retmarks}{Number of marked fish returned to the population}
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Population Size
#'    \item Abundance
#'    \item Mark-Recapture
#'    \item Capture-Recapture
#'    \item Schnabel
#'    \item Schumacher-Eschmeyer
#'  }
#' 
#' @concept Abundance 'Population Size' 'Mark-Recapture' 'Capture-Recapture' Schnabel
#'
#' @source Originally from 
#' 
#' Gerking, S.D.  1953.  Vital statistics of the fish population of Gordy Lake, Indiana.  Transactions of the American Fisheries Society.  82:48-67.
#'
#' But also found in Table 2.4 of 
#'
#' Krebs, C.J.  1999. Ecological Methodology. Addison-Welsey Educational Publishing, second edition.
#'
#' and Table 4.4 of
#'
#' Ricker, W.E.  1975. Computation and interpretation of biological statistics of fish populations. Technical Report Bulletin 191, Bulletin of the Fisheries Research Board of Canada
#'    
#' @keywords datasets
#' 
#' @examples
#' data(SunfishIN)
#' str(SunfishIN)
#' SunfishIN
NULL