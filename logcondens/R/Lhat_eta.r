Lhat_eta <- function(x, w, eta){
    n <- length(x)
    dx <- diff(x)
    phi <- etaphi(x, eta)
    ll <- sum(w * phi) - sum(dx * J00(phi[1:(n - 1)], phi[2:n]))
    return(list(ll = ll))
}
