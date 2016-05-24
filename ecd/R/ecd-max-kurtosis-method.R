#' Utility to calculate where the maximum kurtosis is on the positive j=0 line
#'
#' This utility calculates the kurtosis for alpha from 2.85 to 3.00. Then
#' the location and value of maximum kurtosis is presented.
#' 
#' @param jinv specify 0 (default) or 1728.
#'
#' @return numeric vector, in which the first element is alpha, 
#'         and the second element is the maximum kurtosis.
#'
#' @keywords utility
#'
#' @author Stephen H-T. Lihn
#'
#' @export
#'
#' @examples
#' \dontrun{
#'     k <- ecd.max_kurtosis()
#'     alpha <- k[1]
#'     kurtosis <- k[2]
#' }
### <======================================================================>
"ecd.max_kurtosis" <- function(jinv=0)
{
    if (jinv==0) {
        one = 1 # ecd.mp1
        alpha <- seq(2.85, 3.0, by=0.01)
        K1 <- function(alpha) {
            ecd.kurtosis(ecd(alpha, 0, sigma=one))
        }
        kurt <- simplify2array(parallel::mclapply(alpha,K1))
        
        k_max = max(kurt)
        a_max = alpha[kurt==max(kurt)]
        return(c(a_max, k_max))        
    }
    if (jinv==1728) {
        one = 1 # ecd.mp1
        gamma <- seq(1.480, 1.500, by=0.001)
        K2 <- function(gamma) {
            ecd.kurtosis(ecd(0, gamma, sigma=one))
        }
        kurt <- simplify2array(parallel::mclapply(gamma,K2))
        
        k_max = max(kurt)
        r_max = gamma[kurt==max(kurt)]
        return(c(r_max, k_max))
    }
    stop(paste("Unsupported jinv=", jinv))
}
### <---------------------------------------------------------------------->
