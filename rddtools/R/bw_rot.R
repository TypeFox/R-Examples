#' Bandwidth selector
#' 
#' implements dpill
#' 
#' @param object object of class rdd_data
#' @references McCrary, Justin. (2008) 'Manipulation of the running variable in the regression discontinuity design: A density test,' \emph{Journal of Econometrics}. 142(2): 698-714. \url{http://dx.doi.org/10.1016/j.jeconom.2007.05.005}
#' @export
#' @examples
#' #No discontinuity

### Crary bw

rot_bw <- function(object) {
    
    if (!inherits(object, "rdd_data")) 
        stop("Only works for rdd_data objects")
    cutpoint <- getCutpoint(object)
    x <- object$x
    y <- object$y
    
    ##### first step
    n <- length(y)
    sd_x <- sd(x, na.rm = TRUE)
    bw_pilot <- (2 * sd_x)/sqrt(n)
    
    ## hist
    his <- plotBin(x = x, y = y, h = bw_pilot, cutpoint = cutpoint, plot = FALSE, type = "number")
    # his2 <- hist(x, breaks=c(min(x), his[['x']], max(x)))
    x1 <- his$x
    y1 <- his[, "y.Freq"]
    
    ##### second step
    
    ## regs:
    reg_left <- lm(y1 ~ poly(x1, degree = 4, raw = TRUE), subset = x1 < cutpoint)
    reg_right <- lm(y1 ~ poly(x1, degree = 4, raw = TRUE), subset = x1 >= cutpoint)
    
    
    
}


#' Global bandwidth selector of Ruppert, Sheather and Wand (1995) from package \pkg{KernSmooth}
#' 
#' Uses the global bandwidth selector of Ruppert, Sheather and Wand (1995) 
#' either to the whole function, or to the functions below and above the cutpoint. 
#' 
#' @param object object of class rdd_data created by \code{\link{rdd_data}}
#' @param type Whether to choose a global bandwidth for the whole function (\code{global}) 
#' or for each side (\code{sided})
#' @return One (or two for \code{sided}) bandwidth value. 
#' @references See \code{\link[KernSmooth]{dpill}}
#' @seealso \code{\link{rdd_bw_ik}} Local RDD bandwidth selector using the plug-in method of Imbens and Kalyanaraman (2012)
#' @import KernSmooth
#' @export
#' @examples
#' data(house)
#' rd<- rdd_data(x=house$x, y=house$y, cutpoint=0)
#' rdd_bw_rsw(rd)


rdd_bw_rsw <- function(object, type = c("global", "sided")) {
    
    type <- match.arg(type)
    
    if (!inherits(object, "rdd_data")) 
        stop("Only works for rdd_data objects")
    cutpoint <- getCutpoint(object)
    x <- object$x
    y <- object$y
    
    if (type == "global") {
        bw <- dpill(x = x, y = y)
    } else {
        dat_left <- subset(object, x < cutpoint)
        dat_right <- subset(object, x >= cutpoint)
        
        bw_left <- dpill(x = dat_left$x, y = dat_left$y)
        bw_right <- dpill(x = dat_right$x, y = dat_right$y)
        bw <- c(bw_left, bw_right)
    }
    
    ## result
    bw
} 
