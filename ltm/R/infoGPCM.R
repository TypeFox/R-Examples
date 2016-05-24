infoGPCM <-
function (betas, z, IRT.param) {
    n <- length(z)
    p <- length(betas)
    alphas <- sapply(betas, tail, 1)
    prs <- crf.GPCM(betas, z, IRT.param)
    T.bar <- sapply(prs, function (x) colSums(x * 1:nrow(x)))
    info <- matrix(0, n, p)
    for (j in 1:p) {
        ii <- outer(seq(1, nrow(prs[[j]])), T.bar[, j], "-")^2
        info[, j] <- alphas[j]^2 * colSums(prs[[j]] * ii)
    }
    colnames(info) <- names(betas)
    info
}
