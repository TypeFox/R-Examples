##' A function to return the maximum percentage of the cumulative sum represented 
##' by a subset of the vector
##' @description Returns the proportion of the cumulative sum represented by the 
##' number of elements in the vector a user specifies. This allows the user to 
##' identify the maximum proportion of the total that only X number of elements 
##' may represent in the vector.
##' @param x a numeric vector, missing values are allowed
##' @param cutoff numeric, the number of elements to look at
##' @param na.rm logical, should missing values be excluded?
##' @details Calculates the proportion of a numeric vector reached after sorting 
##' the vector in ascending order and stopping at the specified count
##' @return A numeric proportion
##' @seealso \code{\link{cutoff}} which this function is related to
##' @author Jared E. Knowles
##' @export
##' @examples
##' # for vector
##' a <- rnorm(100, mean=6, sd=1)
##' thresh(a, 8) #return minimum number of elements to account 70 percent of total
thresh <- function(x, cutoff, na.rm = TRUE){
  if(cutoff < 1 | is.na(cutoff)){
    stop("Cutoff value must be an integer 1 or greater and nonmissing")
  }
  if(na.rm){
    x <- x[order(-x)] # sort vector descending
    xb <- cumsum(x)   # take cumulative sum (order matter)
    xc <- xb / sum(x, na.rm=TRUE) # express proportionally
    xc[cutoff]            # report cumulative percentage at given threshold
  } else{
    x <- x[order(-x)] 
    xb <- cumsum(x)   
    xc <- xb / sum(x, na.rm=FALSE)
    xc <- na.omit(xc)
    if(length(xc) == 0){
      NA
    } else {
      xc[cutoff] 
    }
  }
}
