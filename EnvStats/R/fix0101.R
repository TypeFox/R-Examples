fix0101 <-
function (x, alpha = 4, beta = 4, iter = 60, converge = 1e-05) 
{
    if (var(x[, 2]) > 1e-10) {
        cat("warning - censoring values are not all the same", 
            "\n")
    }
    realvalues <- x[x[, 2] < x[, 1], 1]
    pqlvalues <- x[x[, 2] >= x[, 1], 2]
    if (length(realvalues) == 0) 
        cat("warning - no observed x, algorithm will not work", 
            "\n")
    if (length(pqlvalues) == 0) 
        cat("warning - no censored x, algorithm will not work", 
            "\n")
    N <- nrow(x)
    m.little <- length(realvalues)
    n.little <- length(pqlvalues)
    u <- mean(x[, 1])
    s <- sqrt(var(x[, 1]))
    for (k in 1:iter) {
        pql.z <- (pqlvalues - u)/s
        real.z <- (realvalues - u)/s
        gam.wt <- (alpha + 0.5)/(beta + 0.5 * real.z^2)
        real.wt <- sum(realvalues * gam.wt)
        g21 <- g2.m.singly.censored(pql.z, 1, alpha, beta)
        g20 <- g2.m.singly.censored(pql.z, 0, alpha, beta)
        g1.5 <- g1.m.singly.censored(pql.z, 0.5, alpha, beta)
        pql.wt <- u * sum((g21[g20 > 0]/g20[g20 > 0])) - s * 
            sum(g1.5[g20 > 0]/g20[g20 > 0])
        marg.wt <- sum(gam.wt) + sum(g21[g20 > 0]/g20[g20 > 0])
        newmu <- (real.wt + pql.wt)/marg.wt
        real.ss <- sum((realvalues^2 - newmu^2) * gam.wt)
        expected.ss <- (u^2 - newmu^2) * sum(g21[g20 > 0]/g20[g20 > 
            0]) - s * sum(((pqlvalues[g20 > 0] + u) * g1.5[g20 > 
            0])/g20[g20 > 0]) + n.little * s^2
        if (is.na(expected.ss)) {
            cat("warning - fix0101 expected.ss set to zero from NA")
            expected.ss <- 0
        }
        newsig <- sqrt((real.ss + expected.ss)/N)
        eps <- sum(abs(c(u - newmu, s - newsig)))
        iterno <- k
        u <- newmu
        s <- newsig
        if (is.na(eps)) {
            cat("warning fix0101 eps==NA")
            cat("eps = ", eps, " u = ", u, " s = ", s, " real.ss = ", 
                real.ss, " expected.ss = ", expected.ss)
            break
        }
        if (eps < converge) 
            break
    }
    pql.z <- (pqlvalues - u)/s
    real.z <- (realvalues - u)/s
    psi <- sum((real.z * g1.m.singly.censored(real.z, 1.5, alpha, 
        beta))/g1.m.singly.censored(real.z, 0.5, alpha, beta)) - 
        sum((g1.m.singly.censored(pql.z, 0.5, alpha, beta))/g2.m.singly.censored(pql.z, 
            0, alpha, beta))
    chi <- sum(((real.z^2) * g1.m.singly.censored(real.z, 1.5, 
        alpha, beta))/g1.m.singly.censored(real.z, 0.5, alpha, 
        beta)) - sum((pql.z * g1.m.singly.censored(pql.z, 0.5, 
        alpha, beta))/g2.m.singly.censored(pql.z, 0, alpha, beta)) - 
        m.little
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
