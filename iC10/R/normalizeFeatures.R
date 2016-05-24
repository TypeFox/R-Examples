normalizeFeatures <-
function(x, method=c("none", "scale")) {
    if (length(method)==2) method <- "scale"
    if (is.null(x$Exp)) {
        stop("No expression features\n")
    }
    if (method == "none") {
       res <- x
    }
    if (method == "scale") {
        res <- x
        dn <- dimnames(res$Exp)
        res$Exp <- t(apply(res$Exp, 1, scale))
        dimnames(res$Exp) <- dn
    }
#     if (method == "eb") {
#         tmp <- eb(x$train.Exp, x$Exp)
#         res <- x
#         res$train.Exp <- as.matrix(tmp$x)
#         res$Exp <- as.matrix(tmp$y)
#     }
#     if (method == "gq") {
#         tmp <- gq(x$train.Exp, x$Exp)
#         res <- x
#         res$train.Exp <- as.matrix(tmp$x)
#         res$Exp <- as.matrix(tmp$y)
#     }
#     if (method == "xpn") {
#         tmp <- xpn(x$train.Exp, x$Exp)
#         res <- x
#         res$train.Exp <- as.matrix(tmp$x)
#         res$Exp <- as.matrix(tmp$y)
#     }
#     if (method == "train.Metabric") {
#         tmp <- cbind(x$train.Exp, x$Exp)
#         tmp <- normalize.quantiles(tmp)
#         tmp <- t(apply(tmp, 1, scale))
#         tmp1 <- tmp[,1:997]
#         dimnames(tmp1) <- dimnames(x$train.Exp)
#         tmp2 <- tmp[,998:ncol(tmp)]
#         dimnames(tmp2) <- dimnames(x$Exp)
#         res <- x
# 	res$train.Exp <- tmp1
#         res$Exp <- tmp2
#     }
    res
}
