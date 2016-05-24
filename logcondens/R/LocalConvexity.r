LocalConvexity <- function(x, phi){
    n <- length(x)
    deriv <- diff(phi)/diff(x)
    conv <- 1:n * 0
    conv[2:(n - 1)] <- diff(deriv[1:(n - 1)])
    return(matrix(conv, ncol = 1))
}
