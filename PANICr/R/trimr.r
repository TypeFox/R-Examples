
#'@title Trim dataset
#'
#'@description This function trims the dataset by n and/or t rows and columns
#'
#'@usage trimr(x,a,b)
#'
#'@param x A NxT matrix containing the data to be trimmed
#'
#'@param a number of columns to trim from matrix
#'
#'@param b number of rows to trim from matrix
#'
#'@return xx the trimmed data set
#'

trimr <- function(x, a, b) {
    
    x <- as.matrix(x)
    
    nt <- dim(x)[1]
    
    nc <- dim(x)[2]
    
    a0 <- a + 1
    
    b0 <- nt - b
    
    if (a > 0) {
        
        xx <- x[a0:nt, ]
    }
    
    if (b > 0) {
        
        if (a > 0) {
            
            x = xx
        }
        
        xx <- x[1:b0, 1:nc]
    }
    
    return(xx)
    
} 
