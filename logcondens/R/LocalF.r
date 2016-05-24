LocalF <- function(x, phi){
    n <- length(x)
    Fhat <- 1:n * 0
    dx <- diff(x)
    Fhat[2:n] <- cumsum(dx * J00(phi[1:(n - 1)], phi[2:n]))
    return(matrix(Fhat, ncol = 1))
}
