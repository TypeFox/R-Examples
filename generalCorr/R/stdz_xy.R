#' Standardize x and y vectors to force zero mean and unit variance.
#'
#' @param x vector of data which can have NA's
#' @param y  vector of data which can have NA's
#' @importFrom stats sd
#' @return
#' \item{stdx}{standardized values of x}
#' \item{stdy}{standardized values of y}
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @note This works even if there are missing x or y values.
#' @examples
#'
#' \dontrun{
#' set.seed(30)
#' x=sample(20:30)
#' y=sample(21:31)
#' stdz_xy(x,y) }
#'
#' @export

stdz_xy <-
function(x, y){
stdx=(x-mean(x,na.rm=TRUE))/sd(x, na.rm=TRUE)
stdy= (y-mean(y, na.rm=TRUE))/sd(y, na.rm=TRUE)
list(stdx=stdx,stdy=stdy)}
