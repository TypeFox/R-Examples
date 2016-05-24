crf.GPCM2 <-
function (betas, z, IRT.param = TRUE, index, log = FALSE, eps = .Machine$double.eps^(1/2)) {
    eta <- linpred.GPCM(betas, z, IRT.param)
    p <- length(eta)
    out <- vector("list", p)
    for (i in 1:p) {
        num <- exp(rowsum(eta[[i]][index[[i]][, 2], , drop = FALSE], index[[i]][, 1], reorder = FALSE))
        if (!is.matrix(num))
            num <- t(num)
        den <- 1 + colSums(num)
        out[[i]] <- rbind(1/den, num/rep(den, each = nrow(eta[[i]])))
        if (log)
            out[[i]] <- log(out[[i]])
        dimnames(out[[i]]) <- NULL

    }
    out
}
