ics <-
function (icFIT, group, rho, alternative = "two.sided", tol.svd= 10^-8) {
    x <- group
    A <- icFIT$A
    n <- dim(A)[[1]]
    model.int <- 1+(rho>0)
    if (is.numeric(x)) {
        x <- matrix(x, n, 1)
        if (length(unique(x)) == 2)
            test <- "2-sample"
        else test <- "correlation"
    }
    else if (is.character(x) | is.factor(x)) {
        ux <- unique(x)
        nx <- length(ux)
        xout <- matrix(0, n, nx)
        for (i in 1:nx) {
            xout[x == ux[i], i] <- 1
        }
        x <- xout
        test <- paste(nx, "-sample", sep = "")
    }
    q <- dim(x)[[2]]
    pf <- icFIT$pf
    gg <- apply(A, 1, function(x) {
        sum(x * pf)
    })
    m <- length(pf) - 1
    S <- 1 - cumsum(pf[-(m + 1)])
    lambda <- pf[-(m + 1)]/c(1, S[-m])
    Lambda <- cumsum(lambda)
    dpdb.nox <- switch(model.int, S * log(S), -(1/rho)*S * (1 - S^rho))
    calcScores <- function(Ai, mult) {
        sum((Ai[-1] - Ai[-(m + 1)]) * mult)
    }
    dgdb.nox <- apply(A, 1, calcScores, mult = dpdb.nox)
    cc <- dgdb.nox/gg
    mmUpperTri <- matrix(1, m, m)
    mmUpperTri[lower.tri(mmUpperTri)] <- 0
    dpdgam <- switch(model.int, diag(S * log(S)), -diag((1/rho)*S * (1 -
        S^rho)))
    dgdgam <- matrix(NA, n, m)
    for (j in 1:m) {
        dgdgam[, j] <- apply(A, 1, calcScores, mult = dpdgam[j,
            ])
    }
    d2pdb2.nox <- switch(model.int, S * log(S) * (1 + log(S)), -(1/(rho^2))*S * (1 -
        S^rho) * (-1 + (rho+1) * S^rho))
    d2gdb2.nox <- apply(A, 1, calcScores, mult = d2pdb2.nox)
    d2pdbdgam.nox <- switch(model.int, diag(S * log(S) *
        (1 + log(S))), diag(-(1/(rho^2))*S * (1 - S^rho) * (-1 + (rho+1) * S^rho)))
    d2gdbdgam.nox <- matrix(NA, n, m)
    for (j in 1:m) {
        d2gdbdgam.nox[, j] <- apply(A, 1, calcScores, mult = d2pdbdgam.nox[j,
            ])
    }
    d2gdgam2 <- array(0, c(n, m, m))
    if (model.int == 1) {
        for (u in 1:m) {
            d2gdgam2[, u, u] <- (A[, u + 1] - A[, u]) * S[u] *
                log(S[u]) * (1 + log(S[u]))
        }
    }
    else if (model.int == 2) {
        for (u in 1:m) {
            d2gdgam2[, u, u] <- (A[, u + 1] - A[, u]) * (-(1/(rho^2))*S[u]) *
                (1 - (S[u])^rho) * (-1 + (rho+1) * (S[u])^rho)
        }
    }
    else stop("model.int must be 1, or 2")
    d2L.dB2 <- matrix(0, q, q)
    d2L.dgam2 <- matrix(0, m, m)
    d2L.dBdgam <- matrix(0, q, m)
    U <- rep(0, q)
    for (i in 1:n) {
        d2L.dB2 <- d2L.dB2 + (1/gg[i]) * (matrix(x[i, ], q, 1) %*%
            matrix(x[i, ], 1, q)) * (d2gdb2.nox[i] - (1/gg[i]) *
            (dgdb.nox[i]^2))
        d2L.dBdgam <- d2L.dBdgam + (1/gg[i]) * (matrix(x[i, ],
            q, 1) %*% matrix(d2gdbdgam.nox[i, ] - (1/gg[i]) *
            dgdb.nox[i] * dgdgam[i, ], 1, m))
        d2L.dgam2 <- d2L.dgam2 + (1/gg[i]) * (d2gdgam2[i, , ] -
            (1/gg[i]) * matrix(dgdgam[i, ], m, 1) %*% matrix(dgdgam[i,
                ], 1, m))
        U <- U + x[i, ] * cc[i]
    }
    V <- -(d2L.dB2 - d2L.dBdgam %*% solve(d2L.dgam2) %*% t(d2L.dBdgam))
    if (q == 1 | q == 2) {
        Z <- U[1]/sqrt(V[1])
        p.lte <- pnorm(Z)
        p.gte <- 1 - pnorm(Z)
        p.twosidedAbs <- 1 - pchisq(Z^2, 1)
        p.values <- c(p.twosided = p.twosidedAbs, p.lte = p.lte,
            p.gte = p.gte, p.twosidedAbs = p.twosidedAbs)
        if (alternative == "less" | alternative == "greater") {
            statistic <- Z
            names(statistic) <- "Z"
            parameter <- NULL
        }
        else {
            statistic <- Z^2
            names(statistic) <- "Chi Square"
            parameter <- 1
            names(parameter) <- "df"
        }
        p.value <- switch(alternative, less = p.lte, greater = p.gte,
            two.sided = p.twosidedAbs, two.sidedAbs = p.twosidedAbs)
    }
    else {
        svdv <- svd(V)
        index <- (svdv$d > tol.svd)
        ginvV <- svdv$v[, index] %*% ((1/svdv$d[index]) * t(svdv$u[,
            index]))
        chisq.value <- matrix(U, 1, q) %*% ginvV %*% matrix(U,
            q, 1)
        df <- q - 1
        p.twosided <- 1 - pchisq(chisq.value, df)
        p.values <- c(p.twosided = p.twosided, p.twosidedAbs = p.twosided)
        p.value <- p.twosided
        statistic <- chisq.value
        names(statistic) <- "Chi Square"
        parameter <- df
        names(parameter) <- "df"
    }
    out <- list(p.value = p.value, p.values = p.values, statistic = statistic,
        parameter = parameter, V = V, d2L.dB2 = d2L.dB2, d2L.dgam2 = d2L.dgam2,
        d2L.dBdgam = d2L.dBdgam)
    out
}
