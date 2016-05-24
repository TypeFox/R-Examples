#' Given two endpoints gives max and min length of catenary joining them
#'
#' Takes data frame giving endpoints and return min length and max
#' length before hits the ground
#'
#' @param endpoints dataframe with x and y column
#' @return vector of min and max length
#' @author Jono Tuke <simon.tuke@@adelaide.edu.au>
#' @export
#' @note February 11 2013
#' @keywords internal
#' @examples
#' x <- c(-1,1)
#' y <- c(2,2)
#' endpoints <- data.frame(x=x,y=y)
#' minmaxLength(endpoints)
minmaxLength <- function(endpoints){
  if(nrow(endpoints) != 2 | ncol(endpoints) != 2){
    stop("endpoints must be a 2x2 matrix or data frame")
  }
  
  minL <- sqrt(sum((endpoints[1,] - endpoints[2,])^2))
  cat <- catenary(endpoints=endpoints,type='max')
  maxL <- L(cat)
  tab <- c(minL,maxL)
  names(tab) <- c("min","max")
  return(tab)
}
