scoretpm <-
function (thetas, type, constraint, max.guessing) {
    thetas. <- thetas
    thetas <- thetas.tpm(thetas, type, constraint, p)
    betas <- thetas[, 2:3]
    cs <- plogis(thetas[, 1]) * max.guessing
    cs.mat <- matrix(cs, k, p, TRUE)
    pr <- probs(Z %*% t(betas))
    prc <- cs.mat + (1 - cs.mat) * pr
    p.xz <- exp(X %*% t(log(prc)) + mX %*% t(log(1 - prc)))
    p.x <- c(p.xz %*% GHw)
    p.zx <- p.xz / p.x
    pqrc <- prc * (1 - prc)
    pqr <- pr * (1 - pr)
    pqc <- cs * (1 - cs)
    out <- matrix(0, p, 3)
    for (i in 1:p) {
        ypi <- X[, i] - matrix(prc[, i], n, k, TRUE)
        ypi[na.ind[, i], ] <- 0
        mat <- (ypi * p.zx) %*% (GHw * (1 - cs[i]) * pqr[, i] * Z / pqrc[, i])
        if (max.guessing == 1) {
            out[i, 1] <- sum(rep((ypi * p.zx) %*% (GHw * pqc[i] * (1 - pr[, i]) / pqrc[, i]), obs))
        }
        out[i, 2:3] <- colSums(mat[rep(1:n, obs), ])
    }
    if (max.guessing < 1) {
        out[, 1] <- -cd.tpm(thetas., logLiktpm, type = type, constraint = constraint, max.guessing = max.guessing, k = p)
    }
    if (type == "rasch")
        out <- c(out[, 1:2], sum(out[, 3]))
    if (!is.null(constraint)) {
        if (type == "rasch" && any(ind <- constraint[, 2] == 3))
            -out[-c((constraint[!ind, 2] - 1) * p + constraint[!ind, 1], length(out))]
        else
            -out[-((constraint[, 2] - 1) * p + constraint[, 1])]
    } else
        -as.vector(out)
}
