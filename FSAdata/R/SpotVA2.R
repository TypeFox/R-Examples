#' @title Ages (subsample) and lengths (all fish) for Spot.
#' 
#' @description As many as 10 fish per 1-inch total length intervals from the \code{\link[FSA]{SpotVA1}} data frame were obtained for age assignment.  The remaining fish in the file were only measured for length (i.e., the ages were deleted).  This data file can be used to demonstrate the use of age-length keys.
#' 
#' @name SpotVA2
#' 
#' @docType data
#' 
#' @format A data frame of 403 observations on the following 2 variables:
#'  \describe{
#'    \item{tl}{Measured total lengths (in inches)}
#'    \item{age}{Ages assigned from examination of otoliths} 
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Age-Length Key 
#'  }
#' 
#' @concept 'Age-Length Key'
#' 
#' @seealso \code{SpotVA1} in \pkg{FSA}.
#' 
#' @source From Table 1 in Chapter 8 (Spot) of the VMRC Final Report on Finfish Ageing, 2002 by the Center for Quantitative Fisheries Ecology at Old Dominion University.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(SpotVA2)
#' str(SpotVA2)
#' head(SpotVA2)
#' 
#' ## Extract the aged sample
#' spot.aged <- subset(SpotVA2,!is.na(age))
#' str(spot.aged)
#' 
#' ## Extract the length sample
#' spot.length <- subset(SpotVA2,is.na(age))
#' str(spot.length)
#' 
NULL