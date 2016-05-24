Local_LL <- function(x, w, phi){
    n <- length(x)
    dx <- diff(x)
    L <- sum(w * phi) - sum(dx * J00(phi[1:(n - 1)], phi[2:n]))
    return(L)
}
