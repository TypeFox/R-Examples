logLiktpm <-
function (thetas, type, constraint, max.guessing) {
    thetas <- thetas.tpm(thetas, type, constraint, p)
    betas <- thetas[, 2:3]
    cs <- plogis(thetas[, 1]) * max.guessing
    cs <- matrix(cs, k, p, TRUE)
    pr <- cs + (1 - cs) * probs(Z %*% t(betas))
    p.xz <- exp(X %*% t(log(pr)) + mX %*% t(log(1 - pr)))
    p.x <- rep(c(p.xz %*% GHw), obs)
    -sum(log(p.x))
}
