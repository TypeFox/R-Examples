##' Get index by specifying cut points
##'
##' Get index by specifying cut points. Assigns index 1 to lower than the first cut point, 2 to larger than the first cut point but lower than the second cut point, and so on. This function is intended to be used for time dependent covariate in Cox regression models to select covariate according to time point. In this case, cut points is a vector of time points, and index is used to select time-specific covariate.
##' @param variable a numeric value or vector to be converted into index.
##' @param cutpoints a numeric vector specifying cut point values.
##' @param toLowerCategory a logical value specifying which to assign a cut point to lower category or upper category index. The default is assigning to lower category.
##' @param checkOrder a logical value specifying whether to check the order of cut points. Checking the order is skipped by default.
##' @return an index vector.
##' @seealso \code{\link{selectByIndex}}, \code{\link{epifit}}
##' @examples
##' # cutpoint   cp(1)  cp(2)  ...   cp(n-1)  cp(n)
##' # index    1   *   2  *    ... (n-1) *  (n) * (n+1)
##' a <- rnorm(100) * 10
##' b <- getIndex(a, cutpoints=c(-2,-1,0,1,2))
##' @export
getIndex <- function(variable, cutpoints=NULL, toLowerCategory=TRUE, checkOrder=FALSE){

  if (!is.numeric(variable))
    stop("variable should be numeric vector")

  if(length(cutpoints) < 1)
    stop("cutpoints must be specified")

  npoints <- length(cutpoints)
  
  if(checkOrder){
    if(sum((cutpoints[-1] - cutpoints[-npoints]) > 0) < (npoints - 1))
      stop("cutpoints must be ordered")
  }

  if(length(variable) == 1){
    
    result <- numeric(1)
    
    for(i in 1:npoints)
      if(variable < npoints[i])
        return(i);
    
    return(npoints+1)
    
  } else {
    result <- numeric(length(variable))
    cutpoints <- c(-Inf, cutpoints, Inf)
    
    if(toLowerCategory){
      for(i in 1:(npoints+1))
        result <- ifelse(cutpoints[i] < variable & variable <= cutpoints[i+1], i, result)    
    } else {
      for(i in 1:(npoints+1))
        result <- ifelse(cutpoints[i] <= variable & variable < cutpoints[i+1], i, result)    
    }
    
    return(result)
  }
}
