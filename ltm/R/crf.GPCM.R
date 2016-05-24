crf.GPCM <-
function (betas, z, IRT.param = TRUE, log = FALSE, eps = .Machine$double.eps^(1/2)) {
    lapply(linpred.GPCM(betas, z, IRT.param), function (x) {
        num <- exp(apply(x, 2, cumsum))
        if (!is.matrix(num))
            num <- t(num)
        den <- 1 + colSums(num)
        out <- rbind(1/den, num/rep(den, each = nrow(x)))
        if (any(ind <- out == 1))
            out[ind] <- 1 - eps
        if (any(ind <- out == 0))
            out[ind] <- eps
        if (log)
            out <- log(out)
        out
    })
}
