estN <-
function (z, method = "Indep", se = FALSE, nboot = 200){
    METHODS = c("Indep", "HSC", "LSC")
    method = match.arg(method, METHODS)
alpha = 0.05
    z = as.vector(unlist(z))
    FUN = function(z) {
        if (method == "Indep") 
            nhat = estN.Indep(z)
        if (method == "HSC") 
            nhat = estN.HSC(z)
        if (method == "LSC") 
            nhat = estN.LSC(z)
        return(c(est = nhat))
    }
    nhat = FUN(z)
    if (se) {
        z_bs = estN.bootstrap(z, round(nhat), nboot)
        se = sd(apply(z_bs, 2, FUN))
        M = sum(z)
        c = exp(qnorm(1 - alpha / 2) * (log(1 + se^2/(nhat - M)^2))^(0.5))
        cil = M + (nhat - M)/c
        ciu = M + (nhat - M) * c
        out = cbind(est = nhat, se = se, cil = cil, ciu = ciu)
        return(out)
    }
    else return(cbind(nhat))
}
