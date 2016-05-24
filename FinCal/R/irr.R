#' Computing IRR, the internal rate of return
#' 
#' @param cf cash flow,the first cash flow is the initial outlay
#' @seealso \code{\link{pv.uneven}}
#' @seealso \code{\link{npv}}
#' @importFrom stats uniroot
#' @export
#' @examples
#' # irr(cf=c(-5, 1.6, 2.4, 2.8))
irr <- function(cf){
        n <- length(cf)
        subcf <- cf[2:n]
        uniroot(function(r) -1 * pv.uneven(r, subcf) + cf[1], interval=c(1e-10,1e10))$root 
}