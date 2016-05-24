#' Absolute residuals of kernel regression of x on y.
#'
#' This calls the \code{kern} function to implement kernel regression
#' with the option residuals=TRUE and returns afbsolute residuals.
#'
#' The first argument is assumed to be the dependent variable.  If
#' \code{abs_res(x,y)} is used, you are regressing x on y (not the usual y on
#' x)
#'
#' @param x {vector of data on the dependent variable}
#' @param y {vector of data on the regressor}
#' @return absolute values of kernel regression residuals are returned.
#' @note This function is intended for internal use.
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
## @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
### @references %% ~put references to the literature/web site here ~
#' @keywords kern regression residuals
#' @examples
#' \dontrun{
#' set.seed(330)
#' x=sample(20:50)
#' y=sample(20:50)
#' abs_res(x,y)
#' }
#' @export

abs_res <-
function(x, y){
kk1=kern(dep.y=x,reg.x=y,residuals=TRUE)
ares=abs(kk1$resid)
return(ares)}
