MixNRMI1 <-
function (x, probs = c(0.025, 0.5, 0.975), Alpha = 1, Beta = 0, 
    Gama = 0.4, distr.k = 1, distr.p0 = 1, asigma = 0.5, bsigma = 0.5, 
    delta = 3, Delta = 2, Meps = 0.01, Nx = 100, Nit = 500, Pbi = 0.1, 
    epsilon = NULL, printtime = TRUE) 
{
    if (is.null(distr.k)) 
        stop("Argument distr.k is NULL. Should be provided. See help for details.")
    if (is.null(distr.p0)) 
        stop("Argument distr.p0 is NULL. Should be provided. See help for details.")
    tInit <- proc.time()
    n <- length(x)
    y <- x
    for (i in 1:(n/2)) {
        y[i] <- mean(x[1:(n/2)])
    }
    for (i in (n/2):n) {
        y[i] <- mean(x[(n/2):n])
    }
    u <- 1
    sigma <- 1
    if (is.null(epsilon)) 
        epsilon <- sd(x)/4
    xx <- seq(min(x) - epsilon, max(x) + epsilon, length = Nx)
    Fxx <- matrix(NA, nrow = Nx, ncol = Nit)
    fx <- matrix(NA, nrow = n, ncol = Nit)
    R <- seq(Nit)
    S <- seq(Nit)
    U <- seq(Nit)
    Nmt <- seq(Nit)
    mu.p0 = mean(x)
    sigma.p0 = sd(x)
    for (j in seq(Nit)) {
        if (floor(j/100) == ceiling(j/100)) 
            cat("MCMC iteration", j, "of", Nit, "\n")
        tt <- comp1(y)
        ystar <- tt$ystar
        nstar <- tt$nstar
        r <- tt$r
        idx <- tt$idx
        if (Gama != 0) 
            u <- gs3(u, n = n, r = r, alpha = Alpha, beta = Beta, 
                gama = Gama, delta = Delta)
        JiC <- MvInv(eps = Meps, u = u, alpha = Alpha, beta = Beta, 
            gama = Gama, N = 50001)
        Nm <- length(JiC)
        TauiC <- rk(Nm, distr = distr.p0, mu = mu.p0, sigma = sigma.p0)
        ystar <- gs4(ystar, x, idx, distr.k = distr.k, sigma.k = sigma, 
            distr.p0 = distr.p0, mu.p0 = mu.p0, sigma.p0 = sigma.p0)
        Jstar <- rgamma(r, nstar - Gama, Beta + u)
        Tau <- c(TauiC, ystar)
        J <- c(JiC, Jstar)
        tt <- gsHP(ystar, r, distr.p0)
        mu.p0 <- tt$mu.py0
        sigma.p0 <- tt$sigma.py0
        y <- fcondYXA(x, distr = distr.k, Tau = Tau, J = J, sigma = sigma)
        Fxx[, j] <- fcondXA(xx, distr = distr.k, Tau = Tau, J = J, 
            sigma = sigma)
        fx[, j] <- fcondXA(x, distr = distr.k, Tau = Tau, J = J, 
            sigma = sigma)
        sigma <- gs5(sigma, x, y, distr = distr.k, asigma = asigma, 
            bsigma = bsigma, delta = delta)
        R[j] <- r
        S[j] <- sigma
        U[j] <- u
        Nmt[j] <- Nm
    }
    biseq <- seq(floor(Pbi * Nit))
    Fxx <- Fxx[, -biseq]
    qx <- as.data.frame(t(apply(Fxx, 1, quantile, probs = probs)))
    names(qx) <- paste("q", probs, sep = "")
    qx <- cbind(mean = apply(Fxx, 1, mean), qx)
    R <- R[-biseq]
    S <- S[-biseq]
    U <- U[-biseq]
    cpo <- 1/apply(1/fx, 1, mean)
    if (printtime) {
        cat(" >>> Total processing time (sec.):\n")
        print(procTime <- proc.time() - tInit)
    }
    return(list(xx = xx, qx = qx, cpo = cpo, R = R, S = S, U = U, 
        Nm = Nmt, Nx = Nx, Nit = Nit, Pbi = Pbi, procTime = procTime))
}
