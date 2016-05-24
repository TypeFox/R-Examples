gr <-
function (bb) {
    gr.ybi <- - crossprod(Z.ind.i, Z.ind.i %*% bb - yi.eta.yxi) / sigma^2
    Yi <- alpha * (eta.si + rowSums(Zsi * rep(bb, each = GKk)))
    gr.tbi <- numeric(ncz)
    for (k in 1:ncz) {
        gr.tbi[k] <- if (d[i]) {
            alpha * Ztime.i[k] - exp(eta.tw1i) * Pi * sum(wk * exp(eta.wsi + Yi) * alpha * Zsi[, k])
        } else {
            - exp(eta.tw1i) * Pi * sum(wk * exp(eta.wsi + eta.si) * alpha * Zsi[, k])
        }
    }
    gr.bi <- if (diag.D) - bb / D else - solve(D, bb)
    - as.vector(gr.ybi + gr.tbi + gr.bi)
}
