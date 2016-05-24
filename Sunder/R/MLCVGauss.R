MLCVGauss <-
function (gen, D_G, D_E, theta.max = c(10, 100 * max(D_G), 100 * 
    max(D_E), 1, 0.9), theta.min = c(0.01, 0.01 * max(D_G), 0.01 * 
    max(D_E), 0.5, 0.5), ntrain = nrow(gen), nresamp = 1) 
{
    nll.GE <- function(theta, g.transp, m.transp, D_E, D_G) {
        nsite <- ncol(g.transp)
        alpha <- theta[1]
        beta_G <- theta[2]
        beta_E <- theta[3]
        gamma <- theta[4]
        delta <- theta[5]
        base = D_G/beta_G + D_E/beta_E
        Cov <- alpha * ((1 - delta) * exp(-(base)^gamma) + delta * 
            diag(nsite))
        ll <- sum(dmnorm(x = g.transp, mean = m.transp, varcov = Cov, 
            log = TRUE))
        return(-ll)
    }
    nll.G <- function(theta, g.transp, m.transp, D_G) {
        nsite <- ncol(g.transp)
        alpha <- theta[1]
        beta_G <- theta[2]
        gamma <- theta[3]
        delta <- theta[4]
        base = D_G/beta_G
        Cov <- alpha * ((1 - delta) * exp(-(base)^gamma) + delta * 
            diag(nsite))
        ll <- sum(dmnorm(x = g.transp, mean = m.transp, varcov = Cov, 
            log = TRUE))
        return(-ll)
    }
    nll.E <- function(theta, g.transp, m.transp, D_E) {
        nsite <- ncol(g.transp)
        alpha <- theta[1]
        beta_E <- theta[2]
        gamma <- theta[3]
        delta <- theta[4]
        base = D_E/beta_E
        Cov <- alpha * ((1 - delta) * exp(-(base)^gamma) + delta * 
            diag(nsite))
        ll <- sum(dmnorm(x = g.transp, mean = m.transp, varcov = Cov, 
            log = TRUE))
        return(-ll)
    }
    nsite <- nrow(gen)
    nloc <- ncol(gen)
    m <- mean(gen)
    m <- matrix(nrow = nsite, ncol = nloc, data = m, byrow = TRUE)
    m.transp <- t(m)
    g.transp <- t(gen)
    ll.GE <- ll.G <- ll.E <- 0
    for (iresamp in 1:nresamp) {
        train <- sample(1:nsite, ntrain)
        val <- setdiff(1:nsite, train)
        init.par <- runif(n = 5, min = theta.min, max = theta.max)
        res <- optim(par = init.par, fn = nll.GE, gr = NULL, 
            g.transp[, train], m.transp[, train], D_E[train, 
                train], D_G[train, train], method = "L-BFGS-B", 
            lower = theta.min, upper = theta.max)
        theta.ml.GE <- res$par
        res <- optim(par = init.par[-3], fn = nll.G, gr = NULL, 
            g.transp[, train], m.transp[, train], D_G[train, 
                train], method = "L-BFGS-B", lower = theta.min[-3], 
            upper = theta.max[-3])
        theta.ml.G <- res$par
        res <- optim(par = init.par[-2], fn = nll.E, gr = NULL, 
            g.transp[, train], m.transp[, train], D_E[train, 
                train], method = "L-BFGS-B", lower = theta.min[-2], 
            upper = theta.max[-2])
        theta.ml.E <- res$par
        if (ntrain < nsite) {
            base = D_G[val, val]/theta.ml.GE[2] + D_E[val, val]/theta.ml.GE[3]
            Cov <- theta.ml.GE[1] * ((1 - theta.ml.GE[5]) * exp(-(base)^theta.ml.GE[4]) + 
                theta.ml.GE[5] * diag(ntrain))
            ll.GE <- ll.GE + sum(dmnorm(x = g.transp[, val], 
                mean = m.transp[, val], varcov = Cov, log = TRUE))
            base <- D_G[val, val]/theta.ml.G[2]
            Cov <- theta.ml.G[1] * ((1 - theta.ml.G[4]) * exp(-(base)^theta.ml.G[3]) + 
                theta.ml.G[4] * diag(ntrain))
            ll.G <- ll.G + sum(dmnorm(x = g.transp[, val], mean = m.transp[, 
                val], varcov = Cov, log = TRUE))
            base <- D_E[val, val]/theta.ml.E[2]
            Cov <- theta.ml.E[1] * ((1 - theta.ml.E[4]) * exp(-(base)^theta.ml.E[3]) + 
                theta.ml.E[4] * diag(ntrain))
            ll.E <- ll.E + sum(dmnorm(x = g.transp[, val], mean = m.transp[, 
                val], varcov = Cov, log = TRUE))
        }
    }
    if (ntrain == nsite) {
        theta <- list(theta.ml.GE = theta.ml.GE, theta.ml.G = theta.ml.G, 
            theta.ml.E = theta.ml.E)
        res <- list(theta = theta)
    }
    if (ntrain < nsite) {
        mod.lik <- c(ll.GE, ll.G, ll.E)
        names(mod.lik) <- c("G+E", "G", "E")
        res <- list(mod.lik = mod.lik)
    }
    return(res)
}
