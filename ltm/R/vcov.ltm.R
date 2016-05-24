vcov.ltm <-
function (object, robust = FALSE, ...) {
    if (!inherits(object, "ltm"))
        stop("Use only with 'ltm' objects.\n")
    inv.hes <- if (robust) {
        score.vec <- function (betas, X, constraint, GH) {
            p <- nrow(betas)
            q. <- ncol(betas)
            pr <- probs(GH$Z %*% t(betas))
            p.xz <- exp(X %*% t(log(pr)) + (1 - X) %*% t(log(1 - pr)))
            p.x <- c(p.xz %*% GH$GHw)
            p.zx <- p.xz / p.x
            Nt <- GH$GHw * colSums(p.zx)
            scores <- matrix(0, p, q.)
            for (i in 1:p) {
                rit <- GH$GHw * colSums(p.zx * X[, i])
                scores[i, ] <- c(crossprod(rit - pr[, i] * Nt, GH$Z))
            }
            if (!is.null(constraint))
                scores[-((constraint[, 2] - 1) * p + constraint[, 1])]
            else
                c(scores)
        }
        X <- object$X
        if (any(is.na(X)))
            stop("currently the robust estimation of standard errors does not allow for missing values")        
        H <- solve(object$hessian)
        n <- nrow(X)
        nb <- nrow(H)
        S <- lapply(1:n, array, data = 0, dim = c(nb, nb))
        for (m in 1:n) {
            sc <- score.vec(object$coef, X[m, , drop = FALSE], object$constraint, object$GH)
            S[[m]] <- outer(sc, sc)
        }
        S <- matSums(S)
        H %*% S %*% H
    } else
        solve(object$hessian)
    p <- nrow(object$coef)
    nams <- paste(rep(object$ltst$nams, each = p), 1:p, sep = ".")
    if (!is.null(constraint <- object$constraint)) {
        nams <- matrix(nams, p)
        nams <- nams[-((constraint[, 2] - 1) * p + constraint[, 1])]
    }
    dimnames(inv.hes) <- list(nams, nams)
    inv.hes
}
