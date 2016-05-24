.cmicalc <-
function (mi) 
{
    cmi <- gamma(mi/2)/(sqrt(mi/2) * gamma((mi - 1)/2))
    is.na <- is.na(cmi)
    cmi[is.na] <- 1 - 3/(4 * mi[is.na] - 1)
    return(cmi)
}
.con.vcov.UN <-
function (vars, cors) 
{
    dims <- length(vars)
    G <- matrix(1, nrow = dims, ncol = dims)
    G[upper.tri(G)] <- cors
    G[lower.tri(G)] <- t(G)[lower.tri(G)]
    H <- diag(sqrt(vars), nrow = dims, ncol = dims)
    H %*% G %*% H
}
.con.vcov.UN.chol <-
function (vars, covs) 
{
    dims <- length(vars)
    G <- matrix(0, nrow = dims, ncol = dims)
    G[upper.tri(G)] <- covs
    diag(G) <- vars
    crossprod(G)
}
.dnchg <-
function (parms, ai, bi, ci, di, X.fit, random, verbose = FALSE, 
    digits = 4, dnchgcalc, dnchgprec, intCtrl) 
{
    p <- ncol(X.fit)
    k <- length(ai)
    b <- parms[seq_len(p)]
    tau2 <- ifelse(random, exp(parms[p + 1]), 0)
    mu.i <- X.fit %*% cbind(b)
    lli <- rep(NA_real_, k)
    if (!random) {
        for (i in seq_len(k)) {
            lli[i] <- log(.dnchgi(logOR = mu.i[i], ai = ai[i], 
                bi = bi[i], ci = ci[i], di = di[i], random = random, 
                dnchgcalc = dnchgcalc, dnchgprec = dnchgprec))
        }
        if (verbose) 
            cat("ll =", formatC(sum(lli), digits = digits, format = "f"), 
                " ", formatC(b, digits = digits, format = "f"), 
                "\n")
    }
    if (random) {
        for (i in seq_len(k)) {
            res <- try(integrate(.dnchgi, lower = intCtrl$lower, 
                upper = intCtrl$upper, ai = ai[i], bi = bi[i], 
                ci = ci[i], di = di[i], mu.i = mu.i[i], tau2 = tau2, 
                random = random, dnchgcalc = dnchgcalc, dnchgprec = dnchgprec, 
                rel.tol = intCtrl$rel.tol, subdivisions = intCtrl$subdivisions, 
                stop.on.error = FALSE), silent = !verbose)
            if (inherits(res, "try-error")) {
                stop(paste0("Could not integrate over density of non-central hypergeometric distribution in study ", 
                  i, "."))
            }
            else {
                if (res$value > 0) {
                  lli[i] <- log(res$value)
                }
                else {
                  lli[i] <- -Inf
                }
            }
        }
        if (verbose) 
            cat("ll = ", formatC(sum(lli), digits = digits, format = "f"), 
                " ", formatC(tau2, digits = digits, format = "f"), 
                " ", formatC(b, digits = digits, format = "f"), 
                "\n")
    }
    return(-sum(lli))
}
.dnchgi <-
function (logOR, ai, bi, ci, di, mu.i, tau2, random, dnchgcalc, 
    dnchgprec) 
{
    k <- length(logOR)
    dnchgi <- rep(NA_real_, k)
    logOR[logOR < log(1e-12)] <- log(1e-12)
    logOR[logOR > log(1e+12)] <- log(1e+12)
    for (i in seq_len(k)) {
        ORi <- exp(logOR[i])
        if (dnchgcalc == "dnoncenhypergeom") {
            res <- try(.dnoncenhypergeom(x = ai, n1 = ai + bi, 
                n2 = ci + di, m1 = ai + ci, psi = ORi))
        }
        else {
            res <- try(BiasedUrn::dFNCHypergeo(x = ai, m1 = ai + 
                bi, m2 = ci + di, n = ai + ci, odds = ORi, precision = dnchgprec))
        }
        if (inherits(res, "try-error")) {
            stop(paste0("Could not compute density of non-central hypergeometric distribution in study ", 
                i, "."))
        }
        else {
            dnchgi[i] <- res
        }
    }
    if (random) 
        dnchgi <- dnchgi * dnorm(logOR, mu.i, sqrt(tau2))
    return(dnchgi)
}
.dnoncenhypergeom <-
function (x = NA, n1, n2, m1, psi) 
{
    mode.compute <- function(n1, n2, m1, psi, ll, uu) {
        a <- psi - 1
        b <- -((m1 + n1 + 2) * psi + n2 - m1)
        c <- psi * (n1 + 1) * (m1 + 1)
        q <- b + sign(b) * sqrt(b * b - 4 * a * c)
        q <- -q/2
        mode <- trunc(c/q)
        if (uu >= mode && mode >= ll) 
            return(mode)
        else return(trunc(q/a))
    }
    r.function <- function(n1, n2, m1, psi, i) {
        (n1 - i + 1) * (m1 - i + 1)/i/(n2 - m1 + i) * psi
    }
    ll <- max(0, m1 - n2)
    uu <- min(n1, m1)
    if (n1 < 0 | n2 < 0) 
        stop("n1 or n2 negative in dnoncenhypergeom().\n")
    if (m1 < 0 | m1 > (n1 + n2)) 
        stop("m1 out of range in dnoncenhypergeom().\n")
    if (psi <= 0) 
        stop("psi [odds ratio] negative in dnoncenhypergeom().\n")
    if (!is.na(x) & (x < ll | x > uu)) 
        stop("x out of bounds in dnoncenhypergeom().\n")
    if (!is.na(x) & length(x) > 1) 
        stop("x neither missing or scalar in dnoncenhypergeom().\n")
    mode <- mode.compute(n1, n2, m1, psi, ll, uu)
    pi <- array(1, uu - ll + 1)
    shift <- 1 - ll
    if (mode < uu) {
        r1 <- r.function(n1, n2, m1, psi, (mode + 1):uu)
        pi[(mode + 1 + shift):(uu + shift)] <- cumprod(r1)
    }
    if (mode > ll) {
        r1 <- 1/r.function(n1, n2, m1, psi, mode:(ll + 1))
        pi[(mode - 1 + shift):(ll + shift)] <- cumprod(r1)
    }
    pi <- pi/sum(pi)
    if (is.na(x)) {
        return(cbind(ll:uu, pi))
    }
    else {
        return(pi[x + shift])
    }
}
.genperms <-
function (k) 
{
    v <- seq_len(k)
    sub <- function(k, v) {
        if (k == 1L) {
            matrix(v, 1, k)
        }
        else {
            X <- NULL
            for (i in seq_len(k)) {
                X <- rbind(X, cbind(v[i], Recall(k - 1, v[-i])))
            }
            X
        }
    }
    return(sub(k, v[seq_len(k)]))
}
.GENQ.func <-
function (tau2val, P, vi, Q, alpha, k, p, getlower, verbose = FALSE, 
    digits = 4) 
{
    S <- diag(sqrt(vi + tau2val), nrow = k, ncol = k)
    lambda <- Re(eigen(S %*% P %*% S, symmetric = TRUE, only.values = TRUE)$values)
    if (getlower) {
        res <- CompQuadForm::farebrother(Q, lambda[1:(k - p)])$res - 
            alpha
    }
    else {
        res <- (1 - CompQuadForm::farebrother(Q, lambda[1:(k - 
            p)])$res) - alpha
    }
    if (verbose) 
        cat("tau2 =", formatC(tau2val, digits = digits, width = digits + 
            4, format = "f"), " objective =", res, "\n")
    return(res)
}
.genuperms <-
function (x) 
{
    z <- NULL
    sub <- function(x, y) {
        len.x <- length(x)
        if (len.x == 0L) {
            return(y)
        }
        else {
            prev.num <- 0
            for (i in seq_len(len.x)) {
                num <- x[i]
                if (num > prev.num) {
                  prev.num <- num
                  z <- rbind(z, Recall(x[-i], c(y, num)))
                }
            }
            return(z)
        }
    }
    return(sub(x, y = NULL))
}
.invcalc <-
function (X, W, k) 
{
    sWX <- sqrt(W) %*% X
    res.qrs <- qr.solve(sWX, diag(k))
    return(tcrossprod(res.qrs))
}
.is.int.func <-
function (x, eps = 1e-08) 
{
    all(abs(x - 1) < eps)
}
.ll.rma.mv <-
function (par, reml, Y, M, X.fit, k, pX, D.S, Z.G1, Z.G2, Z.H1, 
    Z.H2, sigma2.val, tau2.val, rho.val, gamma2.val, phi.val, 
    sigma2s, tau2s, rhos, gamma2s, phis, withS, withG, withH, 
    struct, g.levels.r, h.levels.r, tol, sparse, cholesky, posdefify, 
    vctransf, verbose, very.verbose, digits, REMLf) 
{
    if (withS) {
        if (vctransf) {
            sigma2 <- ifelse(is.na(sigma2.val), exp(par[1:sigma2s]), 
                sigma2.val)
        }
        else {
            sigma2 <- ifelse(is.na(sigma2.val), par[1:sigma2s], 
                sigma2.val)
            sigma2[sigma2 < 0] <- 0
        }
        for (j in seq_len(sigma2s)) {
            M <- M + sigma2[j] * D.S[[j]]
        }
    }
    if (withG) {
        if (cholesky[1]) {
            tau2 <- par[(sigma2s + 1):(sigma2s + tau2s)]
            rho <- par[(sigma2s + tau2s + 1):(sigma2s + tau2s + 
                rhos)]
        }
        else {
            if (vctransf) {
                tau2 <- ifelse(is.na(tau2.val), exp(par[(sigma2s + 
                  1):(sigma2s + tau2s)]), tau2.val)
                rho <- ifelse(is.na(rho.val), transf.ztor(par[(sigma2s + 
                  tau2s + 1):(sigma2s + tau2s + rhos)]), rho.val)
            }
            else {
                tau2 <- ifelse(is.na(tau2.val), par[(sigma2s + 
                  1):(sigma2s + tau2s)], tau2.val)
                rho <- ifelse(is.na(rho.val), par[(sigma2s + 
                  tau2s + 1):(sigma2s + tau2s + rhos)], rho.val)
                tau2[tau2 < 0] <- 0
                rho[rho > 1] <- 1
                rho[rho < -1] <- -1
            }
        }
        ncol.Z.G1 <- ncol(Z.G1)
        if (struct[1] == "CS") {
            G <- matrix(rho * tau2, nrow = ncol.Z.G1, ncol = ncol.Z.G1)
            diag(G) <- tau2
        }
        if (struct[1] == "HCS") {
            G <- matrix(rho, nrow = ncol.Z.G1, ncol = ncol.Z.G1)
            diag(G) <- 1
            G <- diag(sqrt(tau2), nrow = ncol.Z.G1, ncol = ncol.Z.G1) %*% 
                G %*% diag(sqrt(tau2), nrow = ncol.Z.G1, ncol = ncol.Z.G1)
            diag(G) <- tau2
        }
        if (struct[1] == "UN") {
            if (cholesky[1]) {
                G <- .con.vcov.UN.chol(tau2, rho)
                tau2 <- diag(G)
                rho <- cov2cor(G)[upper.tri(G)]
                tau2[!is.na(tau2.val)] <- tau2.val[!is.na(tau2.val)]
                rho[!is.na(rho.val)] <- rho.val[!is.na(rho.val)]
            }
            G <- .con.vcov.UN(tau2, rho)
            if (posdefify) {
                G <- as.matrix(nearPD(G)$mat)
                tau2 <- diag(G)
                rho <- cov2cor(G)[upper.tri(G)]
            }
        }
        if (struct[1] == "ID" || struct[1] == "DIAG") {
            G <- diag(tau2, nrow = ncol.Z.G1, ncol = ncol.Z.G1)
        }
        if (struct[1] == "UNHO") {
            G <- matrix(1, nrow = ncol.Z.G1, ncol = ncol.Z.G1)
            G[upper.tri(G)] <- rho
            G[lower.tri(G)] <- t(G)[lower.tri(G)]
            G <- diag(sqrt(rep(tau2, ncol.Z.G1)), nrow = ncol.Z.G1, 
                ncol = ncol.Z.G1) %*% G %*% diag(sqrt(rep(tau2, 
                ncol.Z.G1)), nrow = ncol.Z.G1, ncol = ncol.Z.G1)
            if (posdefify) {
                G <- as.matrix(nearPD(G, keepDiag = TRUE)$mat)
                tau2 <- G[1, 1]
                rho <- cov2cor(G)[upper.tri(G)]
            }
        }
        if (struct[1] == "AR") {
            if (ncol.Z.G1 > 1) {
                G <- toeplitz(ARMAacf(ar = rho, lag.max = ncol.Z.G1 - 
                  1))
            }
            else {
                G <- diag(1)
            }
            G <- diag(sqrt(rep(tau2, ncol.Z.G1)), nrow = ncol.Z.G1, 
                ncol = ncol.Z.G1) %*% G %*% diag(sqrt(rep(tau2, 
                ncol.Z.G1)), nrow = ncol.Z.G1, ncol = ncol.Z.G1)
            diag(G) <- tau2
        }
        if (struct[1] == "HAR") {
            if (ncol.Z.G1 > 1) {
                G <- toeplitz(ARMAacf(ar = rho, lag.max = ncol.Z.G1 - 
                  1))
            }
            else {
                G <- diag(1)
            }
            G <- diag(sqrt(tau2), nrow = ncol.Z.G1, ncol = ncol.Z.G1) %*% 
                G %*% diag(sqrt(tau2), nrow = ncol.Z.G1, ncol = ncol.Z.G1)
            diag(G) <- tau2
        }
        if (any(g.levels.r)) {
            G[g.levels.r, ] <- 0
            G[, g.levels.r] <- 0
        }
        if (sparse) 
            G <- Matrix(G, sparse = TRUE)
        M <- M + (Z.G1 %*% G %*% t(Z.G1)) * tcrossprod(Z.G2)
    }
    if (withH) {
        if (cholesky[2]) {
            gamma2 <- par[(sigma2s + tau2s + rhos + 1):(sigma2s + 
                tau2s + rhos + gamma2s)]
            phi <- par[(sigma2s + tau2s + rhos + gamma2s + 1):(sigma2s + 
                tau2s + rhos + gamma2s + phis)]
        }
        else {
            if (vctransf) {
                gamma2 <- ifelse(is.na(gamma2.val), exp(par[(sigma2s + 
                  tau2s + rhos + 1):(sigma2s + tau2s + rhos + 
                  gamma2s)]), gamma2.val)
                phi <- ifelse(is.na(phi.val), transf.ztor(par[(sigma2s + 
                  tau2s + rhos + gamma2s + 1):(sigma2s + tau2s + 
                  rhos + gamma2s + phis)]), phi.val)
            }
            else {
                gamma2 <- ifelse(is.na(gamma2.val), par[(sigma2s + 
                  tau2s + rhos + 1):(sigma2s + tau2s + rhos + 
                  gamma2s)], gamma2.val)
                phi <- ifelse(is.na(phi.val), par[(sigma2s + 
                  tau2s + rhos + gamma2s + 1):(sigma2s + tau2s + 
                  rhos + gamma2s + phis)], phi.val)
                gamma2[gamma2 < 0] <- 0
                phi[phi > 1] <- 1
                phi[phi < -1] <- -1
            }
        }
        ncol.Z.H1 <- ncol(Z.H1)
        if (struct[2] == "CS") {
            H <- matrix(phi * gamma2, nrow = ncol.Z.H1, ncol = ncol.Z.H1)
            diag(H) <- gamma2
        }
        if (struct[2] == "HCS") {
            H <- matrix(phi, nrow = ncol.Z.H1, ncol = ncol.Z.H1)
            diag(H) <- 1
            H <- diag(sqrt(gamma2), nrow = ncol.Z.H1, ncol = ncol.Z.H1) %*% 
                H %*% diag(sqrt(gamma2), nrow = ncol.Z.H1, ncol = ncol.Z.H1)
            diag(H) <- gamma2
        }
        if (struct[2] == "UN") {
            if (cholesky[2]) {
                H <- .con.vcov.UN.chol(gamma2, phi)
                gamma2 <- diag(H)
                phi <- cov2cor(H)[upper.tri(H)]
                gamma2[!is.na(gamma2.val)] <- gamma2.val[!is.na(gamma2.val)]
                phi[!is.na(phi.val)] <- phi.val[!is.na(phi.val)]
            }
            H <- .con.vcov.UN(gamma2, phi)
            if (posdefify) {
                H <- as.matrix(nearPD(H)$mat)
                gamma2 <- diag(H)
                phi <- cov2cor(H)[upper.tri(H)]
            }
        }
        if (struct[2] == "ID" || struct[2] == "DIAG") {
            H <- diag(gamma2, nrow = ncol.Z.H1, ncol = ncol.Z.H1)
        }
        if (struct[2] == "UNHO") {
            H <- matrix(1, nrow = ncol.Z.H1, ncol = ncol.Z.H1)
            H[upper.tri(H)] <- phi
            H[lower.tri(H)] <- t(H)[lower.tri(H)]
            H <- diag(sqrt(rep(gamma2, ncol.Z.H1)), nrow = ncol.Z.H1, 
                ncol = ncol.Z.H1) %*% H %*% diag(sqrt(rep(gamma2, 
                ncol.Z.H1)), nrow = ncol.Z.H1, ncol = ncol.Z.H1)
            if (posdefify) {
                H <- as.matrix(nearPD(H, keepDiag = TRUE)$mat)
                gamma2 <- H[1, 1]
                phi <- cov2cor(H)[upper.tri(H)]
            }
        }
        if (struct[2] == "AR") {
            if (ncol.Z.H1 > 1) {
                H <- toeplitz(ARMAacf(ar = phi, lag.max = ncol.Z.H1 - 
                  1))
            }
            else {
                H <- diag(1)
            }
            H <- diag(sqrt(rep(gamma2, ncol.Z.H1)), nrow = ncol.Z.H1, 
                ncol = ncol.Z.H1) %*% H %*% diag(sqrt(rep(gamma2, 
                ncol.Z.H1)), nrow = ncol.Z.H1, ncol = ncol.Z.H1)
            diag(H) <- gamma2
        }
        if (struct[2] == "HAR") {
            if (ncol.Z.H1 > 1) {
                H <- toeplitz(ARMAacf(ar = phi, lag.max = ncol.Z.H1 - 
                  1))
            }
            else {
                H <- diag(1)
            }
            H <- diag(sqrt(gamma2), nrow = ncol.Z.H1, ncol = ncol.Z.H1) %*% 
                H %*% diag(sqrt(gamma2), nrow = ncol.Z.H1, ncol = ncol.Z.H1)
            diag(H) <- gamma2
        }
        if (any(h.levels.r)) {
            H[h.levels.r, ] <- 0
            H[, h.levels.r] <- 0
        }
        if (sparse) 
            H <- Matrix(H, sparse = TRUE)
        M <- M + (Z.H1 %*% H %*% t(Z.H1)) * tcrossprod(Z.H2)
    }
    if (posdefify) 
        M <- as.matrix(nearPD(M)$mat)
    if (verbose) {
        L <- try(chol(M), silent = !verbose)
    }
    else {
        L <- suppressWarnings(try(chol(M), silent = !verbose))
    }
    if (inherits(L, "try-error")) {
        llval <- -Inf
    }
    else {
        W <- chol2inv(L)
        U <- chol(W)
        sX <- U %*% X.fit
        sY <- U %*% Y
        b <- solve(crossprod(sX), crossprod(sX, sY))
        RSS.f <- sum(as.vector(sY - sX %*% b)^2)
        if (reml) {
            llval <- -1/2 * (k - pX) * log(2 * base::pi) + ifelse(REMLf, 
                1/2 * determinant(crossprod(X.fit), logarithm = TRUE)$modulus, 
                0) - 1/2 * determinant(M, logarithm = TRUE)$modulus - 
                1/2 * determinant(crossprod(X.fit, W) %*% X.fit, 
                  logarithm = TRUE)$modulus - 1/2 * RSS.f
        }
        else {
            llval <- -1/2 * (k) * log(2 * base::pi) - 1/2 * determinant(M, 
                logarithm = TRUE)$modulus - 1/2 * RSS.f
        }
    }
    if ((vctransf && verbose) || (!vctransf && very.verbose)) {
        if (withS) 
            cat("sigma2 =", ifelse(is.na(sigma2), NA, paste(formatC(sigma2, 
                digits = digits, format = "f", flag = " "), " ", 
                sep = "")), "  ", sep = "")
        if (withG) 
            cat("tau2 =", ifelse(is.na(tau2), NA, paste(formatC(tau2, 
                digits = digits, format = "f", flag = " "), " ", 
                sep = "")), "  ", sep = "")
        if (withG) 
            cat("rho =", ifelse(is.na(rho), NA, paste(formatC(rho, 
                digits = digits, format = "f", flag = " "), " ", 
                sep = "")), "  ", sep = "")
        if (withH) 
            cat("gamma2 =", ifelse(is.na(gamma2), NA, paste(formatC(gamma2, 
                digits = digits, format = "f", flag = " "), " ", 
                sep = "")), "  ", sep = "")
        if (withH) 
            cat("phi =", ifelse(is.na(phi), NA, paste(formatC(phi, 
                digits = digits, format = "f", flag = " "), " ", 
                sep = "")), "  ", sep = "")
        cat("  ll = ", ifelse(is.na(llval), NA, formatC(llval, 
            digits = digits, format = "f", flag = " ")), sep = "", 
            "\n")
    }
    return(-1 * llval)
}
.ll.rma.tau2 <-
function (par, yi, vi, X, Z, reml, k, p, verbose, digits, REMLf) 
{
    b.tau2 <- par
    tau2 <- exp(c(Z %*% b.tau2))
    wi <- 1/(vi + tau2)
    W <- diag(wi, nrow = k, ncol = k)
    stXWX <- .invcalc(X = X, W = W, k = k)
    b.fe <- stXWX %*% crossprod(X, W) %*% as.matrix(yi)
    RSS.f <- sum(wi * (yi - X %*% b.fe)^2)
    if (!reml) {
        llval <- -1/2 * (k) * log(2 * base::pi) - 1/2 * sum(log(vi + 
            tau2)) - 1/2 * RSS.f
    }
    else {
        llval <- -1/2 * (k - p) * log(2 * base::pi) + ifelse(REMLf, 
            1/2 * determinant(crossprod(X), logarithm = TRUE)$modulus, 
            0) - 1/2 * sum(log(vi + tau2)) - 1/2 * determinant(crossprod(X, 
            W) %*% X, logarithm = TRUE)$modulus - 1/2 * RSS.f
    }
    if (verbose) {
        cat("ll = ", ifelse(is.na(llval), NA, formatC(llval, 
            digits = digits, format = "f", flag = " ")), " ", 
            sep = "")
        cat("b.tau2 =", ifelse(is.na(b.tau2), NA, paste(formatC(b.tau2, 
            digits = digits, format = "f", flag = " "), " ", 
            sep = "")), "\n", sep = "")
    }
    return(-1 * llval)
}
.modfit <-
function (Y, X, V, k, tol = 1e-07) 
{
    eV <- eigen(V, symmetric = TRUE)
    d <- eV$values
    if (any(d <= tol)) {
        stop("Var-cov matrix is not positive definite.")
    }
    A <- diag(1/sqrt(d), nrow = k, ncol = k) %*% t(eV$vectors)
    Ainv <- eV$vectors %*% diag(1/sqrt(d), nrow = k, ncol = k)
    AX <- A %*% X
    res.qrs <- qr.solve(AX, diag(k))
    vb <- tcrossprod(res.qrs)
    b <- qr.solve(AX, A %*% Y)
    res <- list(b = b, vb = vb)
    return(res)
}
.profile.rma.mv <-
function (val, obj, comp, sigma2.pos, tau2.pos, rho.pos, gamma2.pos, 
    phi.pos, parallel = FALSE, CI = FALSE, objective, verbose = FALSE) 
{
    if (parallel == "snow") 
        library(metafor)
    sigma2.arg <- ifelse(obj$vc.fix$sigma2, obj$sigma2, NA)
    tau2.arg <- ifelse(obj$vc.fix$tau2, obj$tau2, NA)
    rho.arg <- ifelse(obj$vc.fix$rho, obj$rho, NA)
    gamma2.arg <- ifelse(obj$vc.fix$gamma2, obj$gamma2, NA)
    phi.arg <- ifelse(obj$vc.fix$phi, obj$phi, NA)
    if (comp == "sigma2") 
        sigma2.arg[sigma2.pos] <- val
    if (comp == "tau2") 
        tau2.arg[tau2.pos] <- val
    if (comp == "rho") 
        rho.arg[rho.pos] <- val
    if (comp == "gamma2") 
        gamma2.arg[gamma2.pos] <- val
    if (comp == "phi") 
        phi.arg[phi.pos] <- val
    res <- try(suppressWarnings(rma.mv(obj$yi, obj$V, obj$W, 
        mods = obj$X, random = obj$random, struct = obj$struct, 
        intercept = FALSE, method = obj$method, tdist = obj$knha, 
        level = obj$level, R = obj$R, Rscale = obj$Rscale, data = obj$mf.r, 
        sigma2 = sigma2.arg, tau2 = tau2.arg, rho = rho.arg, 
        gamma2 = gamma2.arg, phi = phi.arg, control = obj$control)), 
        silent = TRUE)
    if (!CI) {
        if (inherits(res, "try-error")) {
            list(ll = NA, b = matrix(NA, nrow = nrow(obj$b), 
                ncol = 1), ci.lb = rep(NA, length(obj$ci.lb)), 
                ci.ub = rep(NA, length(obj$ci.ub)))
        }
        else {
            list(ll = logLik(res), b = res$b, ci.lb = res$ci.lb, 
                ci.ub = res$ci.ub)
        }
    }
    else {
        if (inherits(res, "try-error")) {
            if (verbose) 
                cat("vc =", formatC(val, digits = obj$digits, 
                  width = obj$digits + 4, format = "f"), " LRT - objective = NA", 
                  "\n")
            stop()
        }
        else {
            difference <- -2 * (logLik(res) - logLik(obj)) - 
                objective
            if (verbose) 
                cat("vc =", formatC(val, digits = obj$digits, 
                  width = obj$digits + 4, format = "f"), " LRT - objective =", 
                  difference, "\n")
            return(difference)
        }
    }
}
.profile.rma.uni <-
function (val, obj, parallel = FALSE, CI = FALSE, objective, 
    verbose = FALSE) 
{
    if (parallel == "snow") 
        library(metafor)
    res <- try(suppressWarnings(rma.uni(obj$yi, obj$vi, weights = obj$weights, 
        mods = obj$X, method = obj$method, weighted = obj$weighted, 
        intercept = FALSE, knha = obj$knha, level = obj$level, 
        control = obj$control, tau2 = val)), silent = TRUE)
    if (!CI) {
        if (inherits(res, "try-error")) {
            list(ll = NA, b = matrix(NA, nrow = nrow(obj$b), 
                ncol = 1), ci.lb = rep(NA, length(obj$ci.lb)), 
                ci.ub = rep(NA, length(obj$ci.ub)))
        }
        else {
            list(ll = logLik(res), b = res$b, ci.lb = res$ci.lb, 
                ci.ub = res$ci.ub)
        }
    }
    else {
        if (inherits(res, "try-error")) {
            if (verbose) 
                cat("tau2 =", formatC(val, digits = obj$digits, 
                  width = obj$digits + 4, format = "f"), " LRT - objective = NA", 
                  "\n")
            stop()
        }
        else {
            difference <- -2 * (logLik(res) - logLik(obj)) - 
                objective
            if (verbose) 
                cat("tau2 =", formatC(val, digits = obj$digits, 
                  width = obj$digits + 4, format = "f"), " LRT - objective =", 
                  difference, "\n")
            return(difference)
        }
    }
}
.pval <-
function (p, digits = 4, showeq = FALSE, sep = "") 
{
    cutoff <- paste(c(".", rep(0, digits - 1), 1), collapse = "")
    ncutoff <- as.numeric(cutoff)
    ifelse(is.na(p), paste0(ifelse(showeq, "=", ""), sep, NA), 
        ifelse(p >= ncutoff, paste0(ifelse(showeq, "=", ""), 
            sep, formatC(p, digits = digits, format = "f")), 
            paste0("<", sep, cutoff)))
}
.QE.func <-
function (tau2val, Y, vi, X, k, objective, verbose = FALSE, digits = 4) 
{
    wi <- 1/(vi + tau2val)
    W <- diag(wi, nrow = k, ncol = k)
    stXWX <- .invcalc(X = X, W = W, k = k)
    P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
    RSS <- crossprod(Y, P) %*% Y
    if (verbose) 
        cat("tau2 =", formatC(tau2val, digits = digits, width = digits + 
            4, format = "f"), " RSS - objective =", c(RSS - objective), 
            "\n")
    return(RSS - objective)
}
.rtet <-
function (ai, bi, ci, di, maxcor = 0.9999) 
{
    if (!requireNamespace("mvtnorm", quietly = TRUE)) 
        stop("Please install the 'mvtnorm' package to compute this measure.")
    fn <- function(par, ai, bi, ci, di, maxcor, fixcut = FALSE) {
        rho <- par[1]
        cut.row <- par[2]
        cut.col <- par[3]
        if (abs(rho) > maxcor) 
            rho <- sign(rho) * maxcor
        if (fixcut) {
            cut.row <- qnorm((ai + bi)/ni)
            cut.col <- qnorm((ai + ci)/ni)
        }
        R <- matrix(c(1, rho, rho, 1), nrow = 2, ncol = 2)
        p.ai <- mvtnorm::pmvnorm(lower = c(-Inf, -Inf), upper = c(cut.col, 
            cut.row), corr = R)
        p.bi <- mvtnorm::pmvnorm(lower = c(cut.col, -Inf), upper = c(+Inf, 
            cut.row), corr = R)
        p.ci <- mvtnorm::pmvnorm(lower = c(-Inf, cut.row), upper = c(cut.col, 
            +Inf), corr = R)
        p.di <- mvtnorm::pmvnorm(lower = c(cut.col, cut.row), 
            upper = c(+Inf, +Inf), corr = R)
        if (any(p.ai <= 0 || p.bi <= 0 || p.ci <= 0 || p.di <= 
            0)) {
            ll <- -Inf
        }
        else {
            ll <- ai * log(p.ai) + bi * log(p.bi) + ci * log(p.ci) + 
                di * log(p.di)
        }
        return(-ll)
    }
    ni <- ai + bi + ci + di
    if ((ai + bi) == 0L || (ci + di) == 0L || (ai + ci) == 0L || 
        (bi + di) == 0L) 
        return(list(yi = 0, vi = Inf))
    if (bi == 0L && ci == 0L) 
        return(list(yi = 1, vi = 0))
    if (ai == 0L && di == 0L) 
        return(list(yi = -1, vi = 0))
    res <- suppressWarnings(try(optimize(fn, interval = c(-1, 
        1), ai = ai, bi = bi, ci = ci, di = di, maxcor = maxcor, 
        fixcut = TRUE), silent = TRUE))
    if (inherits(res, "try-error")) {
        warning("Could not estimate tetrachoric correlation coefficient.")
        return(list(yi = NA, vi = NA))
    }
    res <- try(optim(par = c(res$minimum, qnorm((ai + bi)/ni), 
        qnorm((ai + ci)/ni)), fn, ai = ai, bi = bi, ci = ci, 
        di = di, maxcor = maxcor, fixcut = FALSE, hessian = TRUE), 
        silent = TRUE)
    if (inherits(res, "try-error")) {
        warning("Could not estimate tetrachoric correlation coefficient.")
        return(list(yi = NA, vi = NA))
    }
    vi <- try(chol2inv(chol(res$hessian))[1, 1], silent = TRUE)
    if (inherits(vi, "try-error")) {
        warning("Could not estimate sampling variance of tetrachoric correlation coefficient.")
        vi <- NA
    }
    yi <- res$par[1]
    if (bi == 0 || ci == 0) 
        yi <- 1
    if (ai == 0 || di == 0) 
        yi <- -1
    return(list(yi = yi, vi = vi, sei = sqrt(vi)))
}
.setlab <-
function (measure, transf.char, atransf.char, gentype) 
{
    if (gentype == 1) 
        lab <- "Observed Outcome"
    if (gentype == 2) 
        lab <- "Overall Estimate"
    if (!is.null(measure)) {
        if (measure == "RR") {
            if (transf.char == "FALSE" && atransf.char == "FALSE") {
                lab <- "Log Relative Risk"
            }
            else {
                lab <- "Transformed Log Relative Risk"
                if (atransf.char == "exp" || atransf.char == 
                  "transf.exp.int") 
                  lab <- "Relative Risk (log scale)"
                if (transf.char == "exp" || transf.char == "transf.exp.int") 
                  lab <- "Relative Risk"
            }
        }
        if (is.element(measure, c("OR", "PETO", "D2OR", "D2ORN", 
            "D2ORL"))) {
            if (transf.char == "FALSE" && atransf.char == "FALSE") {
                lab <- "Log Odds Ratio"
            }
            else {
                lab <- "Transformed Log Odds Ratio"
                if (atransf.char == "exp" || atransf.char == 
                  "transf.exp.int") 
                  lab <- "Odds Ratio (log scale)"
                if (transf.char == "exp" || transf.char == "transf.exp.int") 
                  lab <- "Odds Ratio"
            }
        }
        if (measure == "RD") {
            if (transf.char == "FALSE" && atransf.char == "FALSE") {
                lab <- "Risk Difference"
            }
            else {
                lab <- "Transformed Risk Difference"
            }
        }
        if (measure == "AS") {
            if (transf.char == "FALSE" && atransf.char == "FALSE") {
                lab <- "Arcsine Transformed Risk Difference"
            }
            else {
                lab <- "Transformed Arcsine Transformed Risk Difference"
            }
        }
        if (measure == "PHI") {
            if (transf.char == "FALSE" && atransf.char == "FALSE") {
                lab <- "Phi Coefficient"
            }
            else {
                lab <- "Transformed Phi Coefficient"
            }
        }
        if (measure == "YUQ") {
            if (transf.char == "FALSE" && atransf.char == "FALSE") {
                lab <- "Yule's Q"
            }
            else {
                lab <- "Transformed Yule's Q"
            }
        }
        if (measure == "YUY") {
            if (transf.char == "FALSE" && atransf.char == "FALSE") {
                lab <- "Yule's Y"
            }
            else {
                lab <- "Transformed Yule's Y"
            }
        }
        if (measure == "IRR") {
            if (transf.char == "FALSE" && atransf.char == "FALSE") {
                lab <- "Log Incidence Rate Ratio"
            }
            else {
                lab <- "Transformed Log Incidence Relative Risk"
                if (atransf.char == "exp" || atransf.char == 
                  "transf.exp.int") 
                  lab <- "Incidence Rate Ratio (log scale)"
                if (transf.char == "exp" || transf.char == "transf.exp.int") 
                  lab <- "Incidence Rate Ratio"
            }
        }
        if (measure == "IRD") {
            if (transf.char == "FALSE" && atransf.char == "FALSE") {
                lab <- "Incidence Rate Difference"
            }
            else {
                lab <- "Transformed Incidence Rate Difference"
            }
        }
        if (measure == "IRSD") {
            if (transf.char == "FALSE" && atransf.char == "FALSE") {
                lab <- "Square-Root Transformed Incidence Rate Difference"
            }
            else {
                lab <- "Transformed Square-Root Transformed Incidence Rate Difference"
            }
        }
        if (measure == "MD") {
            if (transf.char == "FALSE" && atransf.char == "FALSE") {
                lab <- "Mean Difference"
            }
            else {
                lab <- "Transformed Mean Difference"
            }
        }
        if (is.element(measure, c("SMD", "SMDH", "PBIT", "OR2D", 
            "OR2DN", "OR2DL"))) {
            if (transf.char == "FALSE" && atransf.char == "FALSE") {
                lab <- "Standardized Mean Difference"
            }
            else {
                lab <- "Transformed Standardized Mean Difference"
            }
        }
        if (measure == "ROM") {
            if (transf.char == "FALSE" && atransf.char == "FALSE") {
                lab <- "Log Ratio of Means"
            }
            else {
                lab <- "Transformed Log Ratio of Means"
                if (atransf.char == "exp" || atransf.char == 
                  "transf.exp.int") 
                  lab <- "Ratio of Means (log scale)"
                if (transf.char == "exp" || transf.char == "transf.exp.int") 
                  lab <- "Ratio of Means"
            }
        }
        if (measure == "RPB") {
            if (transf.char == "FALSE" && atransf.char == "FALSE") {
                lab <- "Point-Biserial Correlation"
            }
            else {
                lab <- "Transformed Point-Biserial Correlation"
            }
        }
        if (is.element(measure, c("COR", "UCOR", "RTET", "RBIS"))) {
            if (transf.char == "FALSE" && atransf.char == "FALSE") {
                lab <- "Correlation Coefficient"
            }
            else {
                lab <- "Transformed Correlation Coefficient"
            }
        }
        if (measure == "ZCOR") {
            if (transf.char == "FALSE" && atransf.char == "FALSE") {
                lab <- "Fisher's z Transformed Correlation Coefficient"
            }
            else {
                lab <- "Transformed Fisher's z Transformed Correlation Coefficient"
                if (atransf.char == "transf.ztor" || atransf.char == 
                  "transf.ztor.int") 
                  lab <- "Correlation Coefficient"
                if (transf.char == "transf.ztor" || transf.char == 
                  "transf.ztor.int") 
                  lab <- "Correlation Coefficient"
            }
        }
        if (measure == "PR") {
            if (transf.char == "FALSE" && atransf.char == "FALSE") {
                lab <- "Proportion"
            }
            else {
                lab <- "Transformed Proportion"
            }
        }
        if (measure == "PLN") {
            if (transf.char == "FALSE" && atransf.char == "FALSE") {
                lab <- "Log Proportion"
            }
            else {
                lab <- "Transformed Log Proportion"
                if (atransf.char == "exp" || atransf.char == 
                  "transf.exp.int") 
                  lab <- "Proportion (log scale)"
                if (transf.char == "exp" || transf.char == "transf.exp.int") 
                  lab <- "Proportion"
            }
        }
        if (measure == "PLO") {
            if (transf.char == "FALSE" && atransf.char == "FALSE") {
                lab <- "Log Odds"
            }
            else {
                lab <- "Transformed Log Odds"
                if (atransf.char == "transf.ilogit" || atransf.char == 
                  "transf.ilogit.int" || atransf.char == "plogis") 
                  lab <- "Proportion (logit scale)"
                if (transf.char == "transf.ilogit" || transf.char == 
                  "transf.ilogit.int" || transf.char == "plogis") 
                  lab <- "Proportion"
                if (atransf.char == "exp" || atransf.char == 
                  "transf.exp.int") 
                  lab <- "Odds (log scale)"
                if (transf.char == "exp" || transf.char == "transf.exp.int") 
                  lab <- "Odds"
            }
        }
        if (measure == "PAS") {
            if (transf.char == "FALSE" && atransf.char == "FALSE") {
                lab <- "Arcsine Transformed Proportion"
            }
            else {
                lab <- "Transformed Arcsine Transformed Proportion"
                if (atransf.char == "transf.iarcsin" || atransf.char == 
                  "transf.iarcsin.int") 
                  lab <- "Proportion (arcsine scale)"
                if (transf.char == "transf.iarcsin" || transf.char == 
                  "transf.iarcsin.int") 
                  lab <- "Proportion"
            }
        }
        if (measure == "PFT") {
            if (transf.char == "FALSE" && atransf.char == "FALSE") {
                lab <- "Double Arcsine Transformed Proportion"
            }
            else {
                lab <- "Transformed Double Arcsine Transformed Proportion"
                if (atransf.char == "transf.ift.hm") 
                  lab <- "Proportion"
                if (transf.char == "transf.ift.hm") 
                  lab <- "Proportion"
            }
        }
        if (measure == "IR") {
            if (transf.char == "FALSE" && atransf.char == "FALSE") {
                lab <- "Incidence Rate"
            }
            else {
                lab <- "Transformed Incidence Rate"
            }
        }
        if (measure == "IRLN") {
            if (transf.char == "FALSE" && atransf.char == "FALSE") {
                lab <- "Log Incidence Rate"
            }
            else {
                lab <- "Transformed Log Incidence Rate"
                if (atransf.char == "exp" || atransf.char == 
                  "transf.exp.int") 
                  lab <- "Incidence Rate (log scale)"
                if (transf.char == "exp" || transf.char == "transf.exp.int") 
                  lab <- "Incidence Rate"
            }
        }
        if (measure == "IRS") {
            if (transf.char == "FALSE" && atransf.char == "FALSE") {
                lab <- "Square-Root Transformed Incidence Rate"
            }
            else {
                lab <- "Transformed Square-Root Transformed Incidence Rate"
                if (atransf.char == "transf.isqrt" || atransf.char == 
                  "transf.isqrt.int") 
                  lab <- "Incidence Rate (square-root scale)"
                if (transf.char == "transf.isqrt" || transf.char == 
                  "transf.isqrt.int") 
                  lab <- "Incidence Rate"
            }
        }
        if (measure == "IRFT") {
            if (transf.char == "FALSE" && atransf.char == "FALSE") {
                lab <- "Freeman-Tukey Transformed Incidence Rate"
            }
            else {
                lab <- "Transformed Freeman-Tukey Transformed Incidence Rate"
            }
        }
        if (measure == "MN") {
            if (transf.char == "FALSE" && atransf.char == "FALSE") {
                lab <- "Mean"
            }
            else {
                lab <- "Transformed Mean"
            }
        }
        if (measure == "MC") {
            if (transf.char == "FALSE" && atransf.char == "FALSE") {
                lab <- "Mean Change"
            }
            else {
                lab <- "Transformed Mean Change"
            }
        }
        if (is.element(measure, c("SMCC", "SMCR", "SMCRH"))) {
            if (transf.char == "FALSE" && atransf.char == "FALSE") {
                lab <- "Standardized Mean Change"
            }
            else {
                lab <- "Transformed Standardized Mean Change"
            }
        }
        if (measure == "ROMC") {
            if (transf.char == "FALSE" && atransf.char == "FALSE") {
                lab <- "Log Ratio of Means"
            }
            else {
                lab <- "Transformed Log Ratio of Means"
                if (atransf.char == "exp" || atransf.char == 
                  "transf.exp.int") 
                  lab <- "Ratio of Means (log scale)"
                if (transf.char == "exp" || transf.char == "transf.exp.int") 
                  lab <- "Ratio of Means"
            }
        }
        if (measure == "ARAW") {
            if (transf.char == "FALSE" && atransf.char == "FALSE") {
                lab <- "Coefficient alpha"
            }
            else {
                lab <- "Transformed Coefficient alpha"
            }
        }
        if (measure == "AHW") {
            if (transf.char == "FALSE" && atransf.char == "FALSE") {
                lab <- "Transformed Coefficient alpha"
            }
            else {
                lab <- "Transformed Coefficient alpha"
                if (atransf.char == "transf.iahw") 
                  lab <- "Coefficient alpha"
                if (transf.char == "transf.iahw") 
                  lab <- "Coefficient alpha"
            }
        }
        if (measure == "ABT") {
            if (transf.char == "FALSE" && atransf.char == "FALSE") {
                lab <- "Transformed Coefficient alpha"
            }
            else {
                lab <- "Transformed Coefficient alpha"
                if (atransf.char == "transf.iabt") 
                  lab <- "Coefficient alpha"
                if (transf.char == "transf.iabt") 
                  lab <- "Coefficient alpha"
            }
        }
    }
    return(lab)
}
.tr <-
function (X) 
return(sum(diag(X)))
