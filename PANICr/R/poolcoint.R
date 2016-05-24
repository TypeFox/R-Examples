#'@title Pooling Function for Cointegration test PANIC (2004)
#'
#'@description This function find the P values for the pooled cointegration test in PANIC (2010)
#'
#'@usage poolcoint(a,x,r)
#'
#'@param a A supplied matrix containg the containing p values
#'
#'@param x The adf test to be pooled
#'
#'@param r The number of factors determined by getnfac()
#'
#'@return pvala a P-value for the pooled cointegration test
#'
#'@return pvalb a P-value for the pooled cointegration test
#'
poolcoint <- function(a, x, r) {
    
    x <- as.matrix(x)
    
    N <- ncol(x)
    
    pval <- matrix(0, N, 1)
    
    r <- ifelse(r > 4, 4, r)
    
    for (i in 1:N) {
        
        aa <- abs(a[, r] - x[, i])
        
        j1 <- min(aa)
        
        j2 <- which.min(aa)
        
        pval[i] <- a[j2, 4]
    }
    
    pvala <- -2 * sum(log(pval))
    
    pvalb <- (pvala - 2 * N)/sqrt(4 * N)
    
    results <- list(pvala = pvala, pvalb = pvalb)
    
    return(results)
    
} 
