#' @title Catches in removal events of Cutthroate Trout and Steelhead of various sizes in two reaches of McGarvey Creek (CA).
#' 
#' @description Catches in removal events of Cutthroate Trout (\emph{Oncorhynchus clarki}) and Steelhead  (\emph{Oncorhynchus mykiss}) of various sizes in two reaches of McGarvey Creek (CA).
#' 
#' @details Sampling was conducted using a Smith Root model 15-D POW electrofisher. Block nets were placed at upstream and downstream reach boundaries. Efforts were made to keep the effort consistent between passes. Fixed electrofisher settings were used to maintain capture probabilities during sampling. Index reaches were rested for at least 90 minutes between passes to allow recovery time for fish not captured.  Fish were measured (fork length in mm) and weighed to the nearest 0.1 gm.  Scales were collected from below the dorsal fin on both left and right sides of selected fish.  After data were recorded for each pass, fish were placed in a loating live car until all sampling was completed.  Fish were then released throughout the reach.
#' 
#' @name SalmonidsMCCA
#' 
#' @docType data
#' 
#' @format A data frame of 5 observations on the following 5 variables:
#'  \describe{
#'    \item{reach}{Sampling location.}
#'    \item{group}{Size or species caught (\code{fry}=both age-0 Cutthroat Trout and Steelhead, \code{Steelhead}=age-1+ Steelhead, or \code{Cutthroat}=age-1+ Cutthroat Trout).}
#'    \item{pass1}{Catch on the first removal pass.}
#'    \item{pass2}{Catch on the second removal pass.}
#'    \item{pass3}{Catch on the third removal pass.}
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Population size
#'    \item Abundance
#'    \item Removal
#'  }
#' 
#' @concept Abundance 'Population Size' Removal
#' 
#' @source From Table 2 of Voight, H.  1999.  Assessment of juvenile salmonid populations in two index reaches of McGarvey Creek, a tributary to the lower Klamath River.  First Year of Investigations - 1998.  [Was (is?) from http://www.krisweb.com/biblio/klamath_yuroktfp_voight_1999_7.pdf.]
#' 
#' @keywords datasets
#' 
#' @examples
#' data(SalmonidsMCCA)
#' str(SalmonidsMCCA)
#' head(SalmonidsMCCA)
#' 
#' ## extract data for one reach and group (e.g., 3rd row)
#' SalmonidsMCCA[3,]
#' 
NULL