fix0201.obs <-
function (x, alpha = 2, beta = 1, iter = 60, converge = 1e-04) 
{
    if (var(x[, 2]) > 1e-10) {
        cat("warning - censoring values are not all the same", 
            "\n")
    }
    realvalues <- x[x[, 2] < x[, 1], 1]
    pqlvalues <- x[x[, 2] >= x[, 1], 2]
    if (length(realvalues) == 0) 
        cat("warning - No observed x, estimation will not work.", 
            "\n")
    if (length(pqlvalues) == 0) 
        cat("note: There is no censored data, estimation will proceed.", 
            "\n")
    if (length(pqlvalues) > 0 & length(realvalues) > 0) 
        cat("warning - cen and uncen present, fix0201.obs will not work", 
            "\n")
    N <- nrow(x)
    m.little <- length(realvalues)
    n.little <- length(pqlvalues)
    u <- mean(x[, 1])
    s <- sqrt(var(x[, 1]))
    if (s == 0) {
        cat("error fix0201.obs - all data are the same value = ", 
            u, "\n")
        break
    }
    for (k in 1:iter) {
        real.z <- (realvalues - u)/s
        gam.wt <- (alpha + 0.5)/(beta + 0.5 * real.z^2)
        real.wt <- sum(realvalues * gam.wt)
        marg.wt <- sum(gam.wt)
        newmu <- real.wt/marg.wt
        real.ss <- sum((realvalues - newmu)^2 * gam.wt)
        newsig <- sqrt((real.ss)/N)
        eps <- sum(abs(c(u - newmu, s - newsig)))
        iterno <- k
        u <- newmu
        s <- newsig
        if (eps < converge) 
            break
    }
    real.z <- (realvalues - u)/s
    psi <- sum((real.z * g1.m.singly.censored(real.z, 1.5, alpha, 
        beta))/g1.m.singly.censored(real.z, 0.5, alpha, beta))
    chi <- sum(((real.z^2) * g1.m.singly.censored(real.z, 1.5, 
        alpha, beta))/g1.m.singly.censored(real.z, 0.5, alpha, 
        beta)) - m.little
    labl <- c("maxiter", "iter", "u estim", "s estim", "psi~0", 
        "chi~0")
    mestimate <- matrix(c(iter, iterno, u, s, psi, chi), 1, 6, 
        dimnames = list(iter, labl[1:6]))
    if (iterno == iter) 
        cat("no convergence ", "psi = ", psi, "chi = ", chi, 
            "\n")
    consistent <- fix0111(d = x[1, 2], alpha = alpha, beta = beta, 
        mu = u, sigma = s, eta = u, kappa = s)
    conestimate <- c(eta = consistent[1], kappa = consistent[3], 
        exp.psi = consistent[2], exp.chi = consistent[4])
    list(mestimate = mestimate, conestimate = conestimate)
}
