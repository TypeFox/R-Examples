#' pairsplot of QQ plots
#' @export
#' @param obj dataframe or matrix
#' @param main - title
#' @examples
#'
#'  tmp = matrix(rep((1:100),times = 4) + rnorm(100*4,0,3),ncol=4)
#' pairsQQ( tmp)
#' @seealso \code{\link{qqplot}} and  \code{\link{pairs}}
pairsQQ = function(obj,main=""){
  pairs(obj, panel = function(x,y){
    r <- qqplot(x , y , plot.it = FALSE)
    points(r,pch=".",cex=2)
    abline(0,1,col=2)
  }
  , lower.panel=NULL
  ,main = main
  )
}
