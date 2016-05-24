fn <-
function (bb) {
    eta.yi <- eta.yxi + rowSums(Z.ind.i * rep(bb, each = ni[i])) 
    Yi <- alpha * (eta.si + rowSums(Zsi * rep(bb, each = GKk)))
    log.p.ybi <- sum(dnorm(yi, eta.yi, sigma, log = TRUE))
    log.p.tbi <- if (d[i]) {
        eta.tw1i + eta.tw2i + alpha * (eta.yxT[i] + sum(Ztime.i * bb)) - exp(eta.tw1i) * Pi * sum(wk * exp(eta.wsi + Yi))
    } else {
        - exp(eta.tw1i) * Pi * sum(wk * exp(eta.wsi + alpha * eta.si))
    }
    log.p.bi <- if (diag.D) - 0.5 * crossprod(bb, bb / D)[1, ] else - 0.5 * crossprod(bb, solve(D, bb))[1, ]
    - (log.p.ybi + log.p.tbi + log.p.bi)
}
