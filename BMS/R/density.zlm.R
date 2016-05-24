density.zlm <-
function (x, reg = NULL, addons = "lesz", std.coefs = FALSE, 
    n = 300, plot = TRUE, hnbsteps = 30, addons.lwd = 1.5, ...) 
{
    addons = gsub("E", "", addons, ignore.case = FALSE)
    addons = gsub("S", "", addons, ignore.case = FALSE)
    addons = gsub("b", "", addons, ignore.case = FALSE)
    addons = gsub("m", "", addons, ignore.case = FALSE)
    addons = gsub("p", "", addons, ignore.case = FALSE)
    N = length(x$residuals)
    K = length(x$coefficients) - 1
    tmo = topmod(1, nmaxregressors = K, bbeta = TRUE, liks = x$marg.lik, 
        ncounts = 1, modelbinaries = matrix(rep(1, K), K, 1), 
        betas = matrix(as.vector(x$coefficients[-1]), K), betas2 = matrix(as.vector(x$coef2moments[-1]), 
            K))
    tokenbma = list(info = list(K = K, N = N), arguments = list(), 
        topmod = tmo, start.pos = integer(0), gprior.info = x$gprior.info, 
        X.data = x$model, reg.names = names(x$coefficients)[-1], 
        bms.call = new("call"))
    class(tokenbma) = "bma"
    return(density.bma(tokenbma, reg = reg, addons = addons, 
        std.coefs = std.coefs, n = n, plot = plot, hnbsteps = hnbsteps, 
        addons.lwd = addons.lwd, ...))
}
