Kci <- function (mod1, mod2, correction = "trans", nsim = 99, ngrid = 200, 
    nrep = 1e+05, r = NULL, simu = "both", spctype = 1) 
{
    verbose <- FALSE
    modnamea <- deparse(substitute(mod1))
    modnameb <- deparse(substitute(mod2))
    if (inherits(mod1, "ppm")) {
        I.ppp <- mod1$Q$data
        #lambdaI <- predict(ppm(mod1$Q$data, mod1$trend), type = "trend",  ngrid = ngrid)
	lambdaI <- predict(mod1, type = "trend",  ngrid = ngrid)
        Isim <- "ppm"
        dataname.a <- mod1$call[[2]]
    }
    if (inherits(mod1, "ecespa.minconfit")) {
        I.ppp <- mod1$data
        lambdaI <- mod1$lambda
        Isim <- "spc"
        dataname.a <- mod1$dataname
    }
    if (inherits(mod2, "ppm")) {
        J.ppp <- mod2$Q$data
        #lambdaJ <- predict(ppm(mod2$Q$data, mod2$trend), type = "trend",  ngrid = ngrid)
	lambdaJ <- predict(mod2, type = "trend",  ngrid = ngrid)
        Jsim <- "ppm"
        dataname.b <- mod2$call[[2]]
    }
    if (inherits(mod2, "ecespa.minconfit")) {
        J.ppp <- mod2$data
        lambdaJ <- mod2$lambda
        Jsim <- "spc"
        dataname.b <- mod2$dataname
    }
    Kia <- Kinhom(I.ppp, lambdaI, correction = correction, r = r)
    mi.r <- Kia$r
    Kia <- Kia[[3]]
    Kib <- Kinhom(J.ppp, lambdaJ, correction = correction, r = mi.r)
    Kib <- Kib[[3]]
    IJ.ppp <- superimpose(a = I.ppp, b = J.ppp)
    Kci.ab.o <- Kcross.inhom(IJ.ppp, i = "a", j = "b", lambdaI, 
        lambdaJ, correction = correction, r = mi.r)[[3]]
    Kci.ba.o <- Kcross.inhom(IJ.ppp, i = "b", j = "a", lambdaJ, 
        lambdaI, correction = correction, r = mi.r)[[3]]
    Kia.s <- NULL
    Kib.s <- NULL
    Kci.ab.s <- NULL
    Kci.ba.s <- NULL
    for (i in 1:nsim) {
        progressreport(i, nsim)
        if (Jsim == "ppm") {
            Jsim.ppp <- rmh(mod2, start = list(x.start = J.ppp), 
                control = list(p = 1, nrep = nrep), verbose = verbose)
        }
        else if (Jsim == "spc") 
            Jsim.ppp <- rIPCP(mod2, type = spctype)
        dentro <- !is.na(lambdaJ[Jsim.ppp, drop = FALSE])
        Jsim.ppp <- Jsim.ppp[dentro]
        IJs.ppp <- superimpose(a = I.ppp, b = Jsim.ppp, W = I.ppp$w)
        IsJ.ppp <- IJs.ppp
        if (simu == "both") {
            if (Isim == "ppm") {
                Isim.ppp <- rmh(mod1, start = list(x.start = I.ppp), 
                  control = list(p = 1, nrep = nrep), verbose = verbose)
            }
            else if (Isim == "spc") 
                Isim.ppp <- rIPCP(mod1, type = spctype)
            dentro <- !is.na(lambdaI[Isim.ppp, drop = FALSE])
            Isim.ppp <- Isim.ppp[dentro]
            IsJ.ppp <- superimpose(a = Isim.ppp, b = J.ppp, W = I.ppp$w)
        }
        Kib.s <- cbind(Kib.s, Kinhom(Jsim.ppp, lambdaJ, correction = correction, 
            r = mi.r, nlarge = Inf)[[3]])
        Kci.ab.s <- cbind(Kci.ab.s, Kcross.inhom(IJs.ppp, i = "a", 
            j = "b", lambdaI, lambdaJ, r = mi.r, correction = correction)[[3]])
        Kci.ba.s <- cbind(Kci.ba.s, Kcross.inhom(IsJ.ppp, i = "b", 
            j = "a", lambdaJ, lambdaI, r = mi.r, correction = correction)[[3]])
        if (simu == "both") {
            Kia.s <- cbind(Kia.s, Kinhom(Isim.ppp, lambdaI, correction = correction, 
                r = mi.r, nlarge = Inf)[[3]])
        }
    }
    result <- list(r = mi.r, kia = Kia, kib = Kib, kci.ab.o = Kci.ab.o, 
        kci.ba.o = Kci.ba.o, kci.ab.s = Kci.ab.s, kci.ba.s = Kci.ba.s, 
        kib.s = Kib.s, kia.s = Kia.s, datanamea = dataname.a, 
        datanameb = dataname.b, modnamea = modnamea, modnameb = modnameb, 
        type = "Kci")
    class(result) <- c("ecespa.kci", class(result))
    return(result)
}
