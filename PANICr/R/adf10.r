#'@title ADF test for PANIC (2010)
#'
#'@description
#' This function performs the ADF tests on the idiosyncratic and common components
#' for PANIC (2010).
#'
#'@usage adf(y,k,p)
#'
#'@param y NxT matrix of data.
#'
#'@param k An integer specifying the maximum lag order for individual
#' ADF regressions. Bai and Ng (2004) suggest 4*(Time)^(.25) rounded
#' to the nearest whole number as the maximum number of lags.
#'
#'@param p A binary selection for 0 or 1. p is the order of the determinisitic
#' function in the regression. 0 is for constant only and 1 is to include a trend.
#'
#'@return tstat A vector of t statistics for each column of the input matrix
#'
adf <- function(y, k, p) {
    
    y <- as.matrix(y)
    
    bigt <- dim(y)[1]
    
    bign <- dim(y)[2]
    
    pval <- matrix(0, 1, bign)
    
    tstat <- matrix(0, 1, bign)
    
    constant <- matrix(1, I(bigt - 1), 1)
    
    trend <- seq(1:I(bigt - 1))
    
    for (i in 1:bign) {
        
        dy <- as.matrix(y[2:bigt, i] - y[1:I(bigt - 1), i])
        
        reg <- as.matrix(y[1:I(bigt - 1), i])
        
        for (j in 1:k) {
            
            reg <- cbind(reg, lagn(dy, j))
            
        }
        
        if (p == 0) {
            
            reg <- cbind(reg, constant)
        }
        
        if (p == 1) {
            
            reg <- cbind(reg, constant, trend)
        }
        
        if (k > 0) {
            
            reg <- trimr(reg, k, 0)
            
            dy <- trimr(dy, k, 0)
        }
        
        alpha <- qr.solve(reg, dy)
        
        e <- dy - reg %*% alpha
        
        sig2 <- t(e) %*% e/nrow(as.matrix(dy))
        
        xx <- solve(crossprod(reg))
        
        tstat[i] <- alpha[1]/sqrt(sig2 * xx[1, 1])
        
        
    }
    
    tstat <- t(as.matrix(tstat))
    output <- tstat
    
    return(output)
    
} 
