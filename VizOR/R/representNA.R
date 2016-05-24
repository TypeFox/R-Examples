##' Convert "NA"s to NAs
##' 
##' A basic utility function converting "NA" strings to R's missing type, NA.
##' 
##' @param x A vector of class 'character', which may contain the string "NA" to
##'   be replaced with NA.
##' @return Returns the vector x, with all "NA" elements recoded as NA.
##' @note It is quite possible that this function should be removed, and that the issue of NAs coded
##' as "NA" ought to be resolved by standardized workflows or SOPs for upstream data management.
##' Conversely, if we retain this function, it should be expanded through an optional argument
##' allowing specification of a set of strings to be recoded as NA; this argument might default to
##' c("NA", "."), for example.
##' @author David C. Norris
##' @keywords manip
##' @examples
##' ## TODO: Provide an example
##' @export representNA
representNA <-
function(x){
  x[as.character(x)=="NA"] <- NA
  x
}

