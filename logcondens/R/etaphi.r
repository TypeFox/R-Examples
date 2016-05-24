etaphi <- function(x, eta){
    n <- length(x)
    dx <- c(0, diff(x))
    phi <- eta[1] + c(0, cumsum(dx[2:n] * eta[2:n]))
    return(phi)
}
