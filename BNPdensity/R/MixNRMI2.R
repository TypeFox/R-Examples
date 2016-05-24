MixNRMI2 <-
function (x, probs = c(0.025, 0.5, 0.975), Alpha = 1, Beta = 0, 
    Gama = 0.4, distr.k = 1, distr.py0 = 1, distr.pz0 = 2, mu.pz0 = 3, 
    sigma.pz0 = sqrt(10), delta = 4, kappa = 2, Delta = 2, Meps = 0.01, 
    Nx = 100, Nit = 500, Pbi = 0.1, epsilon = NULL, printtime = TRUE) 
{
    if (is.null(distr.k)) 
        stop("Argument distr.k is NULL. Should be provided. See help for details.")
    if (is.null(distr.py0)) 
        stop("Argument distr.py0 is NULL. Should be provided. See help for details.")
    tInit <- proc.time()
    n <- length(x)
    y <- x
    for (i in 1:(n/2)) {
        y[i] <- mean(x[1:(n/2)])
    }
    for (i in (n/2):n) {
        y[i] <- mean(x[(n/2):n])
    }
    z <- rep(1, n)
    u <- 1
    if (is.null(epsilon)) 
        epsilon <- sd(x)/4
    xx <- seq(min(x) - epsilon, max(x) + epsilon, length = Nx)
    Fxx <- matrix(NA, nrow = Nx, ncol = Nit)
    fx <- matrix(NA, nrow = n, ncol = Nit)
    R <- seq(Nit)
    U <- seq(Nit)
    Nmt <- seq(Nit)
    mu.py0 = mean(x)
    sigma.py0 = sd(x)
    for (j in seq(Nit)) {
        if (floor(j/100) == ceiling(j/100)) 
            cat("MCMC iteration", j, "of", Nit, "\n")
        tt <- comp2(y, z)
        ystar <- tt$ystar
        zstar <- tt$zstar
        nstar <- tt$nstar
        rstar <- tt$rstar
        idx <- tt$idx
        if (Gama != 0) 
            u <- gs3(u, n = n, r = rstar, alpha = Alpha, beta = Beta, 
                gama = Gama, delta = Delta)
        JiC <- MvInv(eps = Meps, u = u, alpha = Alpha, beta = Beta, 
            gama = Gama, N = 50001)
        Nm <- length(JiC)
        TauyC <- rk(Nm, distr = distr.py0, mu = mu.py0, sigma = sigma.py0)
        TauzC <- rk(Nm, distr = distr.pz0, mu = mu.pz0, sigma = sigma.pz0)
        tt <- gsYZstar(ystar, zstar, nstar, rstar, idx, x, delta, 
            kappa, distr.k = distr.k, distr.py0 = distr.py0, 
            mu.py0 = mu.py0, sigma.py0 = sigma.py0, distr.pz0 = distr.pz0, 
            mu.pz0 = mu.pz0, sigma.pz0 = sigma.pz0)
        ystar <- tt$ystar
        zstar <- tt$zstar
        tt <- gsHP(ystar, rstar, distr.py0)
        mu.py0 <- tt$mu.py0
        sigma.py0 <- tt$sigma.py0
        Jstar <- rgamma(rstar, nstar - Gama, Beta + u)
        Tauy <- c(TauyC, ystar)
        Tauz <- c(TauzC, zstar)
        J <- c(JiC, Jstar)
        tt <- fcondYZXA(x, distr = distr.k, Tauy, Tauz, J)
        y <- tt[, 1]
        z <- tt[, 2]
        Fxx[, j] <- fcondXA2(xx, distr = distr.k, Tauy, Tauz, 
            J)
        fx[, j] <- fcondXA2(x, distr = distr.k, Tauy, Tauz, J)
        R[j] <- rstar
        U[j] <- u
        Nmt[j] <- Nm
    }
    biseq <- seq(floor(Pbi * Nit))
    Fxx <- Fxx[, -biseq]
    qx <- as.data.frame(t(apply(Fxx, 1, quantile, probs = probs)))
    names(qx) <- paste("q", probs, sep = "")
    qx <- cbind(mean = apply(Fxx, 1, mean), qx)
    R <- R[-biseq]
    U <- U[-biseq]
    cpo <- 1/apply(1/fx, 1, mean)
    if (printtime) {
        cat(" >>> Total processing time (sec.):\n")
        print(procTime <- proc.time() - tInit)
    }
    return(list(xx = xx, qx = qx, cpo = cpo, R = R, U = U, Nm = Nmt, 
        Nx = Nx, Nit = Nit, Pbi = Pbi, procTime = procTime))
}
