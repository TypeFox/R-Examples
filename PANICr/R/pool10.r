#'@title Pooling Function for PANIC (2010)
#'
#'@description This function find the P values for the pooled test in PANIC (2010)
#'
#'@usage pool(a,x)
#'
#'
#'@param a A supplied matrix containing p values
#'
#'@param x The adf test to be pooled
#'
#'@return adf31a a P-value for the pooled test
#'
#'@return adf31b a P-value for the pooled test
#'

pool <- function(a, x) {
    
    a <- as.matrix(a)
    
    x <- as.matrix(x)
    
    N <- ncol(x)
    
    pval <- matrix(0, N, 1)
    
    for (i in 1:N) {
        
        aa <- abs(a[, 1] - x[, i])
        
        j1 <- min(aa)
        
        j2 <- which.min(aa)
        
        pval[i, ] <- a[j2, 2]
        
    }
    
    pvala <- -2 * sum(log(pval))
    
    pvalb <- (pvala - 2 * N)/sqrt(4 * N)
    
    output <- list(adf31a = pvala, adf31b = pvalb)
    
    return(output)
} 
