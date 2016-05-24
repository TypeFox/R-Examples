hr <- function(bhr, t, eta.ij, x.i){
    hr <- bhr(t)*exp(eta.ij(t = t, x.i = x.i))
    return(hr)}