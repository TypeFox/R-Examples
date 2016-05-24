#' @title Ages and lengths of Lake Erie Freshwater Drum.
#' 
#' @description Assigned ages (from scales) and measured total lengths for each of 1577 Freshwater Drum (\emph{Aplodinotus grunniens}) from Lake Erie.
#' 
#' @name FWDrumLE1
#' 
#' @docType data
#' 
#' @format A data frame with 1577 observations on the following 2 variables.
#'  \describe{
#'    \item{age}{Assigned ages (from scales).} 
#'    \item{tl}{Measured total lengths (mm).} 
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Growth
#'    \item von Bertalanffy
#'  }
#'  
#' @concept Growth 'von Bertalanffy'
#' 
#' @seealso \code{\link{FWDrumLE2}}.
#' 
#' @source Simulated from Table 3 of Bur, M.T.  1984.  Growth, reproduction, mortality, distribution, and biomass of freshwater drum in Lake Erie.  Journal of Great Lakes Research.  10:48-58.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(FWDrumLE1)
#' str(FWDrumLE1)
#' head(FWDrumLE1)
#' plot(tl~age,data=FWDrumLE1)
#' 
NULL