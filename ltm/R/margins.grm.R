margins.grm <-
function (object, type = c("two-way", "three-way"), rule = 3.5, ...) {
    if (!inherits(object, "grm"))
        stop("Use only with 'grm' objects.\n")
    type <- match.arg(type)
    X <- object$X
    nams <- colnames(X)
    n <- nrow(X)
    betas <- object$coef
    p <- length(betas)
    pr <- iprobs(betas, object$GH$Z)
    GHw <- object$GH$GHw
    X <- data.matrix(object$X)[complete.cases(X), ]
    if (type == "two-way") {
        index <- t(combn(p, 2))
        nindex <- nrow(index)
        margins <- vector("list", nindex)
        for (i in 1:nindex) {
            item1 <- index[i, 1]; p1 <- pr[[item1]]; ncp1 <- ncol(p1)
            item2 <- index[i, 2]; p2 <- pr[[item2]]; ncp2 <- ncol(p2)
            obs <- as.matrix(table(X[, item1], X[, item2]))
            ind <- cbind(rep(1:ncp1, each = ncp2), rep(1:ncp2, ncp1))
            pp <- p1[, ind[, 1]] * p2[, ind[, 2]]
            exp <- obs
            exp[ind] <- n * colSums(GHw * pp)
            resid <- (obs - exp)^2/exp
            margins[[i]] <- list(Obs = obs, Exp = exp, Resid = resid, TotalResid = sum(resid), 
                                    rule = rule * ncp1 * ncp2)
        }
    }
    if (type == "three-way") {
        index <- t(combn(p, 3))
        nindex <- nrow(index)
        margins <- vector("list", nindex)
        for (i in 1:nindex) {
            item1 <- index[i, 1]; p1 <- pr[[item1]]; ncp1 <- ncol(p1)
            item2 <- index[i, 2]; p2 <- pr[[item2]]; ncp2 <- ncol(p2)
            item3 <- index[i, 3]; p3 <- pr[[item3]]; ncp3 <- ncol(p3)
            obs <- as.array(table(X[, item1], X[, item2], X[, item3]))
            ind <- cbind(rep(1:ncp1, each = ncp2), rep(1:ncp2, ncp1))
            ind <- cbind(ind[rep(1:nrow(ind), ncp3), ], rep(1:ncp3, each = nrow(ind)))
            pp <- p1[, ind[, 1]] * p2[, ind[, 2]] * p3[, ind[, 3]]
            exp <- obs
            exp[ind] <- n * colSums(GHw * pp)
            resid <- (obs - exp)^2/exp
            margins[[i]] <- list(Obs = obs, Exp = exp, Resid = resid, TotalResid = sum(resid), 
                                    rule = rule * ncp1 * ncp2 * ncp3)            
        }
    }
    out <- list(margins = margins, type = type, nitems = p, names = nams, call = object$call)
    class(out) <- "margins.grm"
    out
}
