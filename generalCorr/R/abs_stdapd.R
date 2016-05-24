#' Absolute values of gradients (apd's) of kernel regressions of x on y when
#' both x and y are standardized.
#'
#' 1) standardize the data to force mean zero and variance unity, 2) kernel
#' regress x on y, with the option `gradients = TRUE' and finally 3) compute
#' the absolute values of gradients
#'
#' The first argument is assumed to be the dependent variable.  If
#' \code{abs_stdapd(x,y)} is used, you are regressing x on y (not the usual y
#' on x). The regressors can be a matrix with 2 or more columns. The missing values
#' are suitably ignored by the standardization.
#'
#' @param x {vector of data on the dependent variable}
#' @param y {data on the regressors which can be a matrix}
#' @importFrom stats sd
#' @return Absolute values of kernel regression gradients are returned after
#' standardizing the data on both sides so that the magnitudes of amorphous
#' partial derivatives (apd's) are comparable between regression of x on y on
#' the one hand and regression of y on x on the other.
## @note %% ~~further notes~~
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#'
#' @keywords kern regression gradients apd
#' @examples
#' \dontrun{
#' set.seed(330)
#' x=sample(20:50)
#' y=sample(20:50)
#' abs_stdapd(x,y)
#' }
#' @export


abs_stdapd <-
function(x, y){
stdx=function(x)(x-mean(x,na.rm=TRUE))/sd(x, na.rm=TRUE)
#allows y to be a matrix
stx=(x-mean(x,na.rm=TRUE))/sd(x, na.rm=TRUE)
p=NCOL(y)
if (p==1) sty=(y-mean(y,na.rm=TRUE))/sd(y, na.rm=TRUE)
if (p>1) sty=apply(y,2,stdx)
kk1=kern(dep.y=as.vector(stx),reg.x=sty,gradients=TRUE)
agrad=abs(kk1$grad)
return(agrad)}
