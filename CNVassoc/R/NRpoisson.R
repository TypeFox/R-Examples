NRpoisson<-
function (y, X, w, beta, variant, tol, max.iter, verbose)
{
    if (missing(variant))
        variant <- !apply(is.na(beta), 1, any)
    beta <- t(apply(beta, 1, function(x) if (any(is.na(x)))
        rep(x[1], length(x))
    else x))
    J <- NCOL(w)
    n <- NROW(y)
    K <- NROW(beta)
    if (max.iter <= 0) {
        derivatives <- hessianPoisson(beta = beta, y = y, w = w,
            X = X, variant = variant)           ##
        S <- derivatives$S
        H <- derivatives$H
    }
    iter <- 1
    error <- Inf
    while (error > tol & iter < max.iter) {
        derivatives <- hessianPoisson(beta = beta, y = y, w = w,
            X = X, variant = variant)       ##
        S <- derivatives$S
        H <- derivatives$H
        beta.old <- beta
        beta.oldv <- matrix2vector(beta.old, variant)
        betav <- as.vector(beta.oldv - qr.solve(H) %*% S)
        beta <- vector2matrix(betav, variant, J)
        error <- max(abs(beta - beta.old))
        colnames(beta) <- paste("clust", 1:J, sep = "")
        rownames(beta) <- dimnames(X)[2][[1]]
        if (verbose) {
            if (iter == 1)
                cat("---- Newton-Raphson procedure ----\n")
            cat("Iter", iter, "Error", error, "\n")
            print(beta)
            cat("\n")
        }
        iter <- iter + 1
    }
    if (iter > max.iter)
        warning("Maximum iterations reached")
    return(list(beta = beta, variant = variant, score = S, hessian = H))
}

