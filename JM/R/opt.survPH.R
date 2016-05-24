opt.survPH <-
function (thetas) {
    gammas <- thetas[seq_len(ncww)]
    alpha <- thetas[ncww + 1]
    eta.tw <- if (!is.null(WW)) as.vector(WW %*% gammas) else rep(0, n)
    eta.t <- eta.tw + alpha * Y
    eta.s <- alpha * Y2
    exp.eta.s <- exp(eta.s)
    log.lambda0T <- log(lambda0[ind.T0])
    log.lambda0T[is.na(log.lambda0T)] <- 0
    log.hazard <- log.lambda0T + eta.t
    S <- matrix(0, n, k)
    S[unq.indT, ] <- rowsum(lambda0[ind.L1] * exp.eta.s, indT)
    log.survival <- - exp(eta.tw) * S
    log.p.tb <- d * log.hazard + log.survival
    p.bytn <- p.byt * log.p.tb
    -sum(p.bytn %*% wGH, na.rm = TRUE)
}
