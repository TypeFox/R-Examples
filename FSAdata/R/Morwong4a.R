#' @title Ages (subsample) and lengths (all fish) for Morwong from Morwong4.
#' 
#' @description A total of 104 fish dispersed proportionately over 1-cm fork length intervals from the \code{\link{Morwong4}} data frame was obtained for age assignment.  The remaining fish in the file were only measured for length (i.e., the ages were deleted).  This data file can be used to demonstrate the use of age-length keys.
#' 
#' @name Morwong4a
#' 
#' @docType data
#' 
#' @format A data frame with 392 observations on the following 2 variables.
#'  \describe{
#'    \item{fl}{Fork lengths (cm)}
#'    \item{age}{Assigned ages}
#'  }
#' 
#' @section Topic(s):
#'  \itemize{ 
#'    \item Age-Length Key
#'  }
#' 
#' @concept 'Age-Length Key'
#' 
#' @seealso \code{\link{Morwong4}}.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(Morwong4a)
#' str(Morwong4a)
#' head(Morwong4a)
#' 
#' ## extract aged sample
#' m4a.aged <- subset(Morwong4a,!is.na(age))
#' str(m4a.aged)
#' 
#' ## extract length sample
#' m4a.length <- subset(Morwong4a,is.na(age))
#' str(m4a.length)
#' 
NULL