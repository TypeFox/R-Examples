#'@title Estimate the nuisance parameters of the error term
#'
#'@description This function estimates the short-run, long-run, and one sided
#' variance of the error term
#'
#'@usage nuisance(res,k)
#'
#'@param res A matrix consisting of the residuals from a factor model.
#'
#'@param k If fixk is 0, then automatic bandwidth selection is performed.
#' Otherwise, the integer placed here will be the selected bandwidth.
#'
#'@return sig2 The vector of short run variances
#'
#'@return omega2 The vector of long run variances
#'
#'@return half The vector of one-sided variances
#'

nuisance <- function(res, k) {
    
    Tn <- dim(res)[1]
    
    N <- dim(res)[2]
    omega2 <- matrix(0, N, 1)
    
    sig2 <- matrix(0, N, 1)
    
    half <- matrix(0, N, 1)
    
    for (i in 1:N) {
        
        NW <- nw(as.matrix(res[, i]), k)
        
        K <- NW$k
        
        Kw <- NW$w
        
        sig2[i] = crossprod(res[, i])/Tn
        
        omega2[i] <- sig2[i]
        
        half[i] = 0
        
        for (j in 1:K) {
            temp <- t(res[1:I(Tn - j - 1), i]) %*% res[I(j + 1):I(Tn - 1), i]/Tn
            
            omega2[i] <- omega2[i] + 2 * Kw[j] * temp
            
            half[i] <- half[i] + Kw[j] * temp
        }
    }
    
    output = list(sig2 = sig2, omega2 = omega2, half = half)
    return(output)
    
} 
