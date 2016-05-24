#' @title Summarized multiple mark-recapture data for Walleye from four lakes in Northern Minnesota.
#' 
#' @description Summary results of capture histories (number captured, number of recaptured fish, and number of unmarked fish that were marked) for Walleye (\emph{Sander vitreus}) collected from four lakes in Northern Minnesota, USA.
#' 
#' @name WalleyeMN06b
#' 
#' @docType data
#' 
#' @format A data frame with 20 observations on the following 5 variables.
#'  \describe{
#'    \item{lake}{Studied lake (\code{Crooked}, \code{Fourmile}, \code{Island}, \code{Tom})}
#'    \item{date}{Capture date}
#'    \item{catch}{Total fish captured in each sample}
#'    \item{recap}{Marked fish captured in each sample}
#'    \item{retMark}{Marked fish returned to the population}
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
#' @concept Abundance 'Population Size' 'Mark-Recapture' 'Capture-Recapture' 'Schnabel'
#' 
#' @source From appendix 2 in Borkholder, B.D., A.J. Edwards, and C. Olson. 2007.  Spring adult and fall juvenile walleye popluation surveys within the 1854 ceded territory of Minnesota, 2006.  Fond du Lac Division of Resource Management, Technical Report 41.  [Was (is?) from http://www.1854treatyauthority.org/cms/files/REP\%20Fish\%20Walleye\%20Survey\%202006.pdf.]
#' 
#' @keywords datasets
#' 
#' @examples
#' data(WalleyeMN06b)
#' str(WalleyeMN06b)
#' head(WalleyeMN06b)
#' 
NULL