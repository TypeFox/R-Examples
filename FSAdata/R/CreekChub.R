#' @title Ages (subsample) and lengths (all fish) for Creek Chub.
#' 
#' @description Ages (subsample) and lengths (all fish) for Creek Chub (\emph{Semotilus atromaculatus}).
#' 
#' @details As many as 10 fish per 10 mm length interval were sampled for age assignment with scales.
#' 
#' @name CreekChub
#' 
#' @docType data
#' 
#' @format A data frame with 218 observations on the following 2 variables.
#'  \describe{
#'    \item{len}{Total length (mm)}
#'    \item{age}{Assigned ages (yrs; from scales)}
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Age-Length Key
#'    \item Growth
#'  }
#'  
#' @concept 'Age-Length Key' 'Growth'
#' 
#' @note Used in the \href{http://derekogle.com/IFAR/}{Introductory Fisheries Analyses with R} book.
#' 
#' @source Recreated (andom digits were added to construct mm lengths from the cm length bins) from information in Box 15.2 of Quist, M.C., Pegg, M.A., and DeVries, D.R. 2012. Age and growth. In Zale, A.V., Parrish, D.L., and Sutton, T.M., editors. Fisheries Techniques, Third Edition, chapter 15, pages 677-731. American Fisheries Society, Bethesda, MD.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(CreekChub)
#' str(CreekChub)
#' head(CreekChub)
#' xtabs(~age,data=CreekChub)
#' plot(len~age,data=CreekChub)
#' 
NULL