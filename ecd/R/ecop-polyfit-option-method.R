#' Poly fit on option prices
#'
#' The poly fits on logarithm of option prices are performed for each side 
#' of the suggested cusp (specified by \code{k.cusp}).
#' This utiility is used mainly to remove the market data noise for 
#' the calculation of log-slope of option prices.
#'
#' @param k numeric, vector of log-strike
#' @param V numeric, vectors of option prices
#' @param k.cusp length-one numeric, the suggested cusp location
#' @param k.new numeric, vector of log-strike to evaluate the poly fit
#' @param degree.left length-one numeric, specifying the degree of poly fit for the left tail
#' @param degree.right length-one numeric, specifying the degree of poly fit for the right tail
#'
#' @return The state prices from the poly fit
#'
#' @keywords ecop
#'
#' @author Stephen H-T. Lihn
#'
#' @export
#'
#' @importFrom stats lm predict
#' 
### <======================================================================>
"ecop.polyfit_option" <- function(k, V, k.cusp, k.new, degree.left=6, degree.right=6) {

    V <- ecd.mp2f(V)
    
    i_p <- which(k >= k.cusp) 
    i_n <- which(k < k.cusp)

    kp <- k[i_p]
    kn <- k[i_n]
    Vp <- V[i_p]
    Vn <- V[i_n]
    
    fitp <- lm( log(Vp) ~ poly(kp, degree=degree.right))
    fitn <- lm( log(Vn) ~ poly(kn, degree=degree.left))
    
    logV.new <- ifelse(k.new >= k.cusp, 
                    predict(fitp, newdata = data.frame(kp=k.new)),
                    predict(fitn, newdata = data.frame(kn=k.new))
                )
    return(exp(logV.new))
}
### <---------------------------------------------------------------------->
