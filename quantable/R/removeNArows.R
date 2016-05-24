#' Removes rows with more than thresh NA's from matrix
#' @export
#' @return matrix
#' @param obj matrix or dataframe
#' @param thresh - maximum number of NA's / row - if more the row will be removed
#' @examples
#'
#' x = matrix(rnorm(10*10),ncol=10)
#' dim(x)
#' x[3,3] = NA
#' x = removeNArows(x)
#' dim(x)
removeNArows <- function(obj, thresh=0 )
{
  x <- apply(obj,1,function(x){sum(is.na(x))})
  obj <- obj[-which(x>thresh),]
}
