unirootf <- function(exit, u, entry, all.bhr, eta.ij, x.i, pme){
    return(cumallcausehr(entry, exit, all.bhr, eta.ij = eta.ij, x.i = x.i, pme = pme) + log(1 - u))}