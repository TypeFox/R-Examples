phieta <- function(x, phi){
    n <- length(x)
    eta <- c(phi[1], diff(phi)/diff(x))
    return(eta)
}
