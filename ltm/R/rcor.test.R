rcor.test <-
function (mat, p.adjust = FALSE, p.adjust.method = "holm", ...) {
    mat <- data.matrix(mat)
    cor.mat <- cor(mat, ...)
    p <- ncol(mat)
    index <- t(combn(p, 2))
    nindex <- nrow(index)
    pvals <- numeric(nindex)
    for (i in 1:nindex) {
        pvals[i] <- cor.test(mat[, index[i, 1]], mat[, index[i, 2]], ...)$p.value
    }
    if (p.adjust)
        pvals <- p.adjust(pvals, method = p.adjust.method)
    out <- list(cor.mat = cor.mat, p.values = cbind(index, pvals))
    class(out) <- "rcor.test"
    out
}
