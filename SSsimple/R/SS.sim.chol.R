SS.sim.chol <-
function(F, H, Q, R.chol, length.out, beta0 = 0) {
    d <- ncol(H)
    n <- nrow(H)
    tau <- length.out
    params <- .internal.chk.mod.params(F, H, Q, R.chol, P0 = NULL, beta0 = beta0, d, n)
    F <- params$F
    H <- params$H
    Q <- params$Q
    beta0 <- params$beta0
    Beta <- matrix(NA, tau, d)
    Z <- matrix(NA, tau, n)
    Y <- matrix(NA, tau, n)
    Eta <- rmvnorm(tau, rep(0, d), (Q))
	
	xperm <- order(attr(R.chol, "pivot"))
	xordpivchol <- R.chol[ , xperm]
	xordpivchol <- matrix(rnorm(tau * n), nrow = tau) %*% xordpivchol
	xmean <- rep(0, n)
	Epsilon <- sweep(xordpivchol, 2, xmean, "+")
	rm(xordpivchol)
	
	
    for (j in 1:tau) {
        if (j == 1) {
            Beta[j, ] <- F %*% beta0 + Eta[j, ]
        }
        else {
            Beta[j, ] <- F %*% Beta[j - 1, ] + Eta[j, ]
        }
        Y[j, ] <- H %*% Beta[j, ]
        Z[j, ] <- Y[j, ] + Epsilon[j, ]
    }
    return(list("Beta" = Beta, "Y" = Y, "Z" = Z))
}

