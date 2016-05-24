`est.sc` <-
function (Yin, fm.X, region, model = "Poi", gmat, nmat, totalit, fm.ga = TRUE, t.i = NULL, phi0 = 1, omega0 = 0, r0 = 1, beta0 = NULL, gamma0 = NULL, sigma0 = 1, psi0 = 1, Tau = 10, alpha = 2){
    if (fm.X == ~1) 
        Xin <- matrix(1, length(Yin), 1)
    else Xin <- model.matrix(fm.X)
    Xin <- cbind(region, Xin)
    if (is.null(beta0) == TRUE) 
        beta0 <- rep(0, (dim(Xin)[2] - 1))
    if (is.null(gamma0) == TRUE) 
        gamma0 <- rep(0, length(gmat))
    if (is.null(t.i) == TRUE) 
        t.i <- rep(1, dim(Yin)[1])
    if (model == "Poi") {
        poi <- poiind(Yin, Xin, t.i, fm.ga, gmat, nmat, totalit, phi0, beta0, gamma0, sigma0, psi0, Tau, alpha)
        return(poi)
    }
    if (model == "NB") {
        nb <- nbind(Yin, Xin, t.i, fm.ga, gmat, nmat, totalit, r0, beta0, gamma0, sigma0, psi0, Tau, alpha)
        return(nb)
    }
    if (model == "GP") {
        gp <- gpind(Yin, Xin, t.i, fm.ga, gmat, nmat, totalit, phi0, beta0, gamma0, sigma0, psi0, Tau, alpha)
        return(gp)
    }
    if (model == "ZIP") {
        zip <- zipind(Yin, Xin, t.i, fm.ga, gmat, nmat, totalit, omega0, beta0, gamma0, sigma0, psi0, Tau, alpha)
        return(zip)
    }
    if (model == "ZIGP") {
        zigp <- zigpind(Yin, Xin, t.i, fm.ga, gmat, nmat, totalit, phi0, omega0, beta0, gamma0, sigma0, psi0, Tau, alpha)
        return(zigp)
    }
}

