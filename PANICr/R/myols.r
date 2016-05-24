#'@title Beta Coefficients for standard OLS
#'
#'@description Returns the Beta values for an Ordinary Least Squares.
#'
#'@usage myols(x,y)
#'
#'@param x A NxT matrix of dependent variables.
#'
#'@param y A 1xT vector that is used as the independent variable.
#'
myols <- function(x, y) {
    
    bhat <- qr.solve(crossprod(x)) %*% t(x) %*% y
    return(bhat)
} 
