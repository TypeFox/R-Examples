loglikltm <-
function (betas, constraint) {
    betas <- betas.ltm(betas, constraint, p, q.)
    pr <- probs(Z %*% t(betas))
    p.xz <- exp(X %*% t(log(pr)) + mX %*% t(log(1-pr)))
    p.x <- rep(c(p.xz %*% GHw), obs)
    -sum(log(p.x))
}
