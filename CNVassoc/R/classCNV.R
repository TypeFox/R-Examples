classCNV <-
function (x, mixture, threshold.CNV.0, threshold.CNV.k)
{
    mm <- mixture$parameter$mu
    ss <- mixture$parameter$sigma
    pp <- mixture$parameter$pi
    k <- length(ss)
    ans <- matrix(NA, nrow = length(x), ncol = k)
    for (i in 1:k) ans[, i] <- pp[i] * dnorm(x, mm[i], ss[i])
    if (any(x <= threshold.CNV.0)) {
        k <- k + 1
        ans <- cbind(ifelse(x <= threshold.CNV.0, Inf, 0), ans)
    }
    if (any(x >= threshold.CNV.k)) {
        k <- k + 1
        ans <- cbind(ans, ifelse(x >= threshold.CNV.k, Inf, 0))
    }
    out <- apply(ans, 1, function(x) which.max(x)[1])
    probs <- ans/rowSums(ans)
    probs[is.na(probs)] <- 1
    attr(out, "probabilities") <- probs
    attr(out, "k") <- k
    out
}
