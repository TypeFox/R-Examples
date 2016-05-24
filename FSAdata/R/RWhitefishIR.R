#' @title Ages and lengths of Round Whitefish.
#'
#' @description Ages and total lengths of Round Whitefish (\emph{Prosopium cylindraecum}) collected from Isle Royale, Lake Superior.
#'
#' @name RWhitefishIR
#' 
#' @docType data
#' 
#' @format A data frame with 103 observations on the following 2 variables.
#'   \describe{
#'     \item{tl}{Total length (mm).}
#'     \item{age}{Age (from scales).} 
#'   }
#'   
#' @section Topic(s):
#'   \itemize{
#'     \item Growth
#'     \item von Bertalanffy
#'   }
#'   
#' @concept Growth 'von Bertalanffy'
#'
#' @source From Table 4 (a random tenths digit was added to the TL) in Bailey, M.M.  1963.  Age, growth, and maturity of round whitefish of the Apostle Islands and Isle Royale Regions, Lake Superior.  Fishery Bulletin, 63:63-75.   [Was (is?) from http://fishbull.noaa.gov/63-1/bailey.pdf.]
#'
#' @keywords datasets
#'
#' @examples
#' data(RWhitefishIR)
#' str(RWhitefishIR)
#' head(RWhitefishIR)
#' plot(tl~age,data=RWhitefishIR)
#'
NULL