`rcusp` <-
function (n, alpha, beta) # with rejection sampling
{
    p = zapsmall(polyroot(c(alpha, beta, 0, -1)))
    q = zapsmall(polyroot(c(beta, alpha^2, 2 * alpha * beta, 
        beta^2 - 3, -2 * alpha, -2 * beta, 0, 1)))
    r = unique(Re(sort(c(p, q))))
    C = dcusp(0, alpha, beta)
    hdcusp = exp((alpha + (beta/2 - r^2/4) * r) * r) * C
    hproposal = hdcusp[2:length(r) - 1 + (hdcusp[-1] > hdcusp[-length(r)])]
    gamma = max(abs(c(alpha, beta)))
    plow = hdcusp[1] * 0.5 * sqrt(pi/gamma)
    phigh = hdcusp[length(r)] * 0.5 * sqrt(pi/gamma)
    pprop = c(plow, hproposal * diff(r), phigh)
    normc = sum(pprop)
    pprop = pprop/normc
    hproposal = hproposal/normc
    r0 = c(min(r), r, max(r))
    rproposal = function(n) {
        nobs = rmultinom(1, n, pprop)
        mr = min(r)
        yl = -abs(rnorm(nobs[1], sd = 1/sqrt(2 * gamma))) + mr
        zl = exp(-gamma * (yl - mr)^2) * hdcusp[1]/normc
        xr = max(r)
        yh = abs(rnorm(nobs[length(nobs)], sd = 1/sqrt(2 * gamma))) + 
            max(r)
        zh = exp(-gamma * (yh - xr)^2) * hdcusp[length(r)]/normc
        y = c()
        z = c()
        for (i in 2:length(r)) {
            y = c(y, runif(nobs[i], r[i - 1], r[i]))
            z = c(z, rep(hproposal[i - 1], nobs[i]))
        }
        list(y = c(yl, y, yh), dens = c(zl, z, zh))
    }
    y = c()
    acc.rate = 0
    while (length(y) < n) {
        .theta = rproposal(n - length(y))
        theta = .theta$y
        .theta$accept = runif(n - length(y)) < exp((alpha + (0.5 * 
            beta - 0.25 * theta^2) * theta) * theta) * C/(normc * 
            .theta$dens)
        if (length(y) < 1) 
            acc.rate = mean(.theta$accept)
        y = c(y, .theta$y[.theta$accept])
        if (length(y) > n) 
            browser()
    }
    y = y[sample(n)]
    y
}

