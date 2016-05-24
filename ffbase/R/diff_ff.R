#' Lagged Differences
#' 
#' Returned suitably lagged and iterated differences
#' 
#' @method diff ff
#' @param x a \code{ff} vector containing values to be differenced 
#' @param lag a n integer indicating which lag to use
#' @param differences an integer indicating the order of the difference
#' @param ... other parameters will be passed on to diff
#' @export
#' @export diff.ff
diff.ff <- function(x, lag=1L, differences = 1L, ...){
  
  d <- NULL  
  
  i.last <- NULL
  for (i in chunk(x, ...)){
    Log$chunk(i)
    i.x <- x[i]
    d <- ffappend(d, diff(c(i.last, i.x), lag=lag))
    i.last <- tail(i.x, lag)
  }
  
  if (differences > 1){
    diff.ff(d, lag=lag, differences=differences-1 , ...)  
  } else {
    d
  }
}


#x<- ff(1:10)
#diff(x, 2)
