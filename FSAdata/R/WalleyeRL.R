#' @title Growth increment data for Red Lakes Walleye.
#' 
#' @description Growth increment data for Red Lakes Walleye (\emph{Sander vitreus}) in one-fish-per-line format.
#' 
#' @note Data is in one-fish-per-line format.
#' 
#' @name WalleyeRL
#' 
#' @docType data
#' 
#' @format A data frame with 1543 observations on the following 13 variables.
#'  \describe{
#'    \item{fish}{A fish identification number.  Unique within a year but not across years.} 
#'    \item{yearcap}{Year the fish was captured.}
#'    \item{ce}{A factor denoting capture gear (\code{C}=commercial and \code{E}=experimental nets).} 
#'    \item{agecap}{Age of fish at capture.}
#'    \item{lencap}{Length of fish at capture.} 
#'    \item{inc1}{Scale measurement to first annulus.} 
#'    \item{inc2}{Scale measurement between first and second annulus.} 
#'    \item{inc3}{Scale measurement between second and third annulus.}
#'    \item{inc4}{Scale measurement between third and fourth annulus.}
#'    \item{inc5}{Scale measurement between fourth and fifth annulus.}
#'    \item{inc6}{Scale measurement between fifth and sixth annulus.}
#'    \item{inc7}{Scale measurement between sixth and seventh annulus.}
#'    \item{radcap}{Scale radius at time of capture} 
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Growth increment analysis 
#'    \item Weisberg linear growth model 
#'  }
#'  
#' @concept Growth Weisberg LGM
#' 
#' @source Cyterski, M.J. and G.R. Spangler.  1996.  A tool for age determination. North American Journal of Fisheries Management, 16:403-412.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(WalleyeRL)
#' str(WalleyeRL)
#' head(WalleyeRL)
#' 
NULL