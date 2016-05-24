select.model2<-
function (pp, sigmas, r, nlarge = 10000, q = 1/4, p = 2, correction = "trans") 
{
    len.sig <- length(sigmas)
    pprho <- pp$n/area.owin(pp$window)
    lower = c(sigma2 = r[2]/10, rho = 1/area.owin(pp$window))
    upper = c(sigma2 = (max(r) * 4)^2, rho = pprho)
    sigma2.0 <- (upper["sigma2"] - lower["sigma2"])/2
    rho.0 <- (upper["rho"] - lower["rho"])/2
    parscale <- c(1, 1)
    HPPs = list()
    models <- list()
    for (i in 1:len.sig) {
        progressreport(i, len.sig)
        lambda <- density.ppp(pp, sigma = sigmas[i], at = "points")
        hpc.model <- ipc.estK2(pp, lambda = lambda, correction = correction, 
            r = r, nlarge = nlarge, p = p, q = q, sigma2 = sigma2.0, 
            rho = rho.0, method = "L-BFGS-B", lower = lower, 
            upper = upper, control = list(parscale = parscale))
        models[[i]] <- lambda
        models[[i + len.sig]] <- hpc.model
        HPPs[[i]] <- Kinhom(pp, lambda = lambda, r = r, correction = correction, 
            nlarge = nlarge)
    }
    pc.model <- ipc.estK2(pp, correction = correction, r = r, 
        nlarge = nlarge, p = p, q = q, sigma2 = sigma2.0, rho = rho.0, 
        method = "L-BFGS-B", lower = lower, upper = upper, control = list(parscale = c(max(upper), 
            min(lower))))
    models[[(2 * len.sig) + 1]] <- pc.model
    homo.lam <- pp$n/area.owin(pp$window)
    P <- Kest(pp, r = r, correction = correction, nlarge = nlarge)
    models[[(2 * len.sig) + 2]] <- homo.lam
    dtheta.fun <- function(Kobs, Kfit) return(sum((Kobs^q - Kfit^q)^p))
    dtheta <- NULL
    Kas <- list()
    for (i in 1:length(sigmas)) {
        dtheta[i] <- dtheta.fun(HPPs[[i]][[3]], pi * HPPs[[i]]$r^2)
        Kas[[i]] <- HPPs[[i]][[3]]
    }
    for (i in 1:length(sigmas)) {
        dtheta[length(sigmas) + i] <- dtheta.fun(models[[length(sigmas) + 
            i]]$Kobs, models[[length(sigmas) + i]]$Kfit)
        Kas[[length(sigmas) + i]] <- models[[length(sigmas) + 
            i]]$Kobs
    }
    dtheta[2 * length(sigmas) + 1] <- dtheta.fun(pc.model$Kobs, 
        pc.model$Kfit)
    dtheta[2 * length(sigmas) + 2] <- dtheta.fun(P[[3]], pi * 
        P$r^2)
    Kas[[2 * length(sigmas) + 1]] <- pc.model$Kobs
    Kas[[2 * length(sigmas) + 2]] <- P[[3]]
    nombres.modelos <- c(paste("HPP_sg_", sigmas), paste("HPC_sg_", 
        sigmas), "PC", "P")
    names(dtheta) <- nombres.modelos
    best.dtheta <- which.min(dtheta)
    names(models) <- nombres.modelos
    nHPP <- length(sigmas)
    nHPC <- nHPP
    nP <- 1
    nPC <- 1
    npar <- c(rep(2, nHPP), rep(4, nHPC), 3, 1)
    aics <- apply(cbind(dtheta, npar), 1, function(x) aic.function(r, 
        x[1], x[2]))
    result <- list(dtheta = dtheta, best.dtheta = dtheta[best.dtheta], 
        best.model = models[[best.dtheta]], models = models, 
        HPPs = HPPs, sigmas = sigmas, aics = aics, Kas = Kas,pp=pp) # incluye el pp original para poder hacer envelopes
    class(result) <- c("selectedmod", class(result))
    return(result)
}