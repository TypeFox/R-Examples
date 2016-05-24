LocalNormalize <- function(x, phi){
    n <- length(x)
    dx <- diff(x)
    phi_new <- phi - log(sum(dx * J00(phi[1:(n - 1)], phi[2:n])))
    return(matrix(phi_new, ncol = 1))
}
