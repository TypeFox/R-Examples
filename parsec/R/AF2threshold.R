AF2threshold <-
function(mpi, prof, zeta = NULL) {
    if (class(mpi)!="ophi") stop("mpi is not an object of class 'ophi'")
    if (class(prof)!="wprof") stop("prof is not an object of class 'wprof'")
    if (!is.null(zeta)) {
        if (class(zeta)!="incidence")
            stop("zeta is not an object of class 'incidence'")
    } else zeta <- getzeta(prof)
    ds <- apply(prof$profiles, 1, mpi$rho_k)
    gen.downset(zeta, ds)
}
