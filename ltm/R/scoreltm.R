scoreltm <-
function (betas, constraint) {
    betas <- betas.ltm(betas, constraint, p, q.)
    pr <- probs(Z %*% t(betas))
    p.xz <- exp(X %*% t(log(pr)) + mX %*% t(log(1 - pr)))
    p.x <- c(p.xz %*% GHw)
    p.zx <- (p.xz / p.x) * obs
    scores <- matrix(0, p, q.)
    for (i in 1:p) {
        ind. <- na.ind[, i]
        Y <- outer(X[, i], pr[, i], "-")
        Y[ind., ] <- 0
        scores[i, ] <- -colSums((p.zx * Y) %*% (Z * GHw))
    }
    if (!is.null(constraint))
        scores[-((constraint[, 2] - 1) * p + constraint[, 1])]
    else
        c(scores)
}
