#' Helper function to obtain running means for prodlim objects.
#' 
#' Compute average values of a variable according to neighborhoods.
#' 
#' 
#' @param x Object of class \code{"neighborhood"}.
#' @param y Vector of numeric values.
#' @param \dots Not used.
#' @author Thomas Alexander Gerds <tag@@biostat.ku.dk>
#' @seealso \code{\link{neighborhood}}
#' @keywords survival
#' @examples
#' 
#' meanNeighbors(x=1:10,y=c(1,10,100,1000,1001,1001,1001,1002,1002,1002))
#'
#' @export
meanNeighbors <- function(x,y,...){
  nnn=neighbors(x,y,...)
  out <- data.frame(x=nnn$nbh$values,
                    y=sapply(nnn$list,mean))
  names(out) <- c("uniqueX","averageY")
  out
}
