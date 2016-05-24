#' Cast for classes derived from survData
#' 
#' Use this function to access \code{survData} methods on an object of a 
#' derived class (e.g. \code{reproData})
#' 
#' @param x an S3 object of a class derived from \code{survData}
#' 
#' @examples 
#' data(zinc)
#' x <- reproData(zinc)
#' 
#' # Compare
#' plot(x)
#'
#' #and
#' plot(as.survData(x))
#' 
#' @export
as.survData <- function(x) {
  if(inherits(x,"survData")) {
    class(x) <- c("survData","data.frame")
    return(x)
  }
  else stop("The class of x is not derived from survData.")
}
