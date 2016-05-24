#' @title Catch-at-age for Bull Trout in Trestle Creek, ID.
#' 
#' @description Catch-at-age (actually carcasses-at-age) for Bull Trout (\emph{Salvelinus confluentis}) in Trestle Creek, ID.
#' 
#' @name BullTroutTC
#' 
#' @docType data
#' 
#' @format A data frame with 6 observations on the following 2 variables.
#'  \describe{
#'    \item{age}{A numeric vector of assigned ages (from otoliths).}
#'    \item{carcasses}{A numeric vector of number of carcasses found in and along Trestle Creek.} 
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Mortality
#'    \item Catch curve
#'  }
#'  
#' @concept Mortality 'Catch Curve'
#' 
#' @source From (approximately) Figure 4a in Downs, C.C., D. Horan, E. Morgan-Harris, and R. Jakubowski.  2006.  Spawning demographics and juvenile dispersal of an adfluvial bull trout population in Trestle Creek, Idaho.  North American Journal of Fisheries Management 26:190-200.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(BullTroutTC)
#' str(BullTroutTC)
#' head(BullTroutTC)
#' plot(log(carcasses)~age,data=BullTroutTC)
#' 
NULL
