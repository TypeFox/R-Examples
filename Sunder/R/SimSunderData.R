SimSunderData <-
function (mod, theta, nsite, nloc, hap.pop.size, nalM, nalm, 
    var.par, scale.par) 
{
    alpha <- theta[1]
    beta_G <- theta[2]
    beta_E <- theta[3]
    gamma <- theta[4]
    delta <- theta[5]
    s <- matrix(nrow = nsite, ncol = 2, runif(2 * nsite))
    D_G <- as.matrix(dist(s))
    Cov.e <- sqrt(var.par) * exp(-D_G/scale.par)
    L = t(chol(Cov.e))
    e <- L %*% matrix(ncol = 1, rnorm(nsite))
    e <- matrix(e)
    D_E <- as.matrix(dist(e))
    if (mod == "G+E") 
        base = D_G/beta_G + D_E/beta_E
    else if (mod == "G") 
        base = D_G/beta_G
    else if (mod == "E") 
        base = D_E/beta_E
    Cov <- (1 - delta) * exp(-(base)^gamma) + delta * diag(nsite)
    theta = c(alpha, beta_G, beta_E, gamma, delta)
    pop = matrix(nrow = nsite, ncol = nloc, data = hap.pop.size)
    if (nalM == 2) 
        nal <- rep(2, nloc)
    else nal <- sample(x = nalm:nalM, size = nloc, replace = TRUE)
    mdim = c(nsite, nloc, nalM)
    f = x = y = gen = array(dim = mdim)
    L = t(chol(Cov))
    z <- array(dim = mdim, rnorm(n = prod(mdim)))
    for (a in 1:nalM) {
        y[, , a] = L %*% z[, , a]
    }
    x = qgamma(pnorm(y), shape = alpha)
    for (iloc in 1:nloc) {
        x[, iloc, -(1:nal[iloc])] = 0
    }
    ss <- matrix(nrow = nsite, ncol = nloc)
    for (i in 1:nsite) {
        X = x[i, , ]
        ss[i, ] = apply(X, 1, sum)
        f[i, , ] = x[i, , ]/ss[i, ]
    }
    for (isite in 1:nsite) {
        for (iloc in 1:nloc) {
            gen[isite, iloc, ] = rmultinom(n = 1, size = pop[isite, 
                iloc], prob = f[isite, iloc, ])
        }
    }
    output <- list(gen = gen, e = e, s = s, D_G = D_G, D_E = D_E, 
        mod = mod, theta = theta)
    return(output)
}
