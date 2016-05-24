#'@title Create lags for matrix
#'
#'@description This function adds lags to a vector.
#'
#'@usage lagn(x,n)
#'
#'@param x A 1XT vector to apply the lag to
#'
#'@param n the number of lags to add to the vector
#'
#'@return y A 1XT vector with a lag n added
#'


lagn <- function(x, n) {
    
    x <- as.matrix(x)
    
    nt <- dim(x)[1]
    
    nc <- dim(x)[2]
    
    if (n > 0) {
        
        x1 <- as.matrix(trimr(x, 0, n))
        
        y <- rbind(matrix(0, n, nc), x1)
    }
    
    if (n < 0) {
        
        x1 <- trimr(x, abs(n), 0)
        
        y <- rbind(x1, matrix(0, abs(n), nc))
    }
    
    
    return(y)
} 
