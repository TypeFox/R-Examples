Score.phGH <-
function (thetas) {
    betas <- thetas[1:ncx]
    sigma <- exp(thetas[ncx + 1])
    gammas <- thetas[seq(ncx + 2, ncx + 1 + ncww)]
    alpha <- thetas[ncx + ncww + 2]
    D <- thetas[seq(ncx + ncww + 3, length(thetas))]
    D <- if (diag.D) exp(D) else chol.transf(D)
    eta.yx <- as.vector(X %*% betas)
    eta.yxT <- as.vector(Xtime %*% betas)
    eta.yxT2 <- as.vector(Xtime2 %*% betas)
    Y <- eta.yxT + Ztime.b
    Y2 <- eta.yxT2 + Ztime2.b
    eta.tw <- if (!is.null(WW)) as.vector(WW %*% gammas) else rep(0, n)
    eta.t <- eta.tw + alpha * Y
    eta.s <- alpha * Y2
    exp.eta.s <- exp(eta.s)
    exp.eta.tw <- exp(eta.tw)
    mu.y <- eta.yx + Ztb
    logNorm <- dnorm(y, mu.y, sigma, TRUE)
    log.p.yb <- rowsum(logNorm, id); dimnames(log.p.yb) <- NULL
    log.lambda0T <- log(lambda0[ind.T0])
    log.lambda0T[is.na(log.lambda0T)] <- 0
    log.hazard <- log.lambda0T + eta.t
    Int <- lambda0[ind.L1] * exp.eta.s
    S <- matrix(0, n, k)
    S[unq.indT, ] <- rowsum(Int, indT, reorder = FALSE)
    log.survival <- - exp.eta.tw * S
    log.p.tb <- d * log.hazard + log.survival
    log.p.b <- if (control$typeGH == "simple") {
        rep(dmvnorm(b, rep(0, ncz), D, TRUE), each = n)
    } else {
        matrix(dmvnorm(do.call(rbind, lis.b), rep(0, ncz), D, TRUE), n, k, byrow = TRUE)
    }
    p.ytb <- exp(log.p.yb + log.p.tb + log.p.b)
    if (control$typeGH != "simple")
        p.ytb <- p.ytb * VCdets
    dimnames(p.ytb) <- NULL
    p.yt <- c(p.ytb %*% wGH)
    p.byt <- p.ytb / p.yt
    post.b <- if (control$typeGH == "simple") {
        p.byt %*% (b * wGH)
    } else {
        sapply(seq_len(ncz), function (i)
            (p.byt * t(sapply(lis.b, "[", seq_len(k), i))) %*% wGH)
    }
    post.vb <- if (control$typeGH == "simple") { 
        if (ncz == 1) {
            c(p.byt %*% (b2 * wGH)) - c(post.b * post.b)
        } else {
            (p.byt %*% (b2 * wGH)) - t(apply(post.b, 1, function (x) x %o% x))
        }
    } else {
        dd <- sapply(seq_len(ncz^2), function (i)
            (p.byt * t(sapply(lis.b2, "[", seq_len(k), i))) %*% wGH)
        bb <- apply(post.b, 1, function (x) x %o% x)
        dd - if (ncz == 1) c(bb) else t(bb)
    }
    Zb <- if (ncz == 1) post.b[id] else rowSums(Z * post.b[id, ], na.rm = TRUE)
    mu <- y - eta.yx
    tr.tZZvarb <- sum(ZtZ * post.vb, na.rm = TRUE)
    sc1 <- - crossprod(X, y - eta.yx - Zb) / sigma^2
    Int.alpha <- Int * alpha
    sc2 <- numeric(ncx)
    for (i in 1:ncx) {
        S1 <- matrix(0, n, k)
        S1[unq.indT, ] <- rowsum(Int.alpha * Xtime2[, i], indT, reorder = FALSE)
        ki <- exp.eta.tw * S1
        kii <- c((p.byt * ki) %*% wGH)
        sc2[i] <- - sum(d * alpha * Xtime[, i] - kii, na.rm = TRUE)
    }
    score.y <- c(c(sc1 + sc2), 
        - sigma * (- N / sigma + drop(crossprod(mu, mu - 2 * Zb) + crossprod(Zb) + tr.tZZvarb) / sigma^3))
    sc.gammas <- if (!is.null(WW))
        - colSums(WW * (d - c((p.byt * (exp.eta.tw * S)) %*% wGH)), na.rm = TRUE)
    else 
        NULL
    S1[unq.indT, ] <- rowsum(Int * Y2, indT, reorder = FALSE)
    sc.alpha <- - sum((p.byt * (d * Y - exp.eta.tw * S1)) %*% wGH, na.rm = TRUE)
    score.t <- c(sc.gammas, sc.alpha)
    score.b <- if (diag.D) {
        svD <- 1 / D
        svD2 <- svD^2
        cS.postVB <- colSums(as.matrix(post.vb), na.rm = TRUE)
        dim(cS.postVB) <- c(ncz, ncz)
        D * 0.5 * (n * svD - diag(cS.postVB) * svD2 - colSums(as.matrix(post.b^2), na.rm = TRUE) * svD2)
    } else {
        svD <- solve(D)
        dD <- deriv.D(D)
        ndD <- length(dD)
        D1 <- sapply(dD, function (x) sum(svD * x))
        D2 <- t(sapply(dD, function (x) c(svD %*% x %*% svD)))
        cS.postVB <- colSums(as.matrix(post.vb), na.rm = TRUE)
        out <- numeric(ndD)
        for (i in seq_along(dD)) {
            D.mat <- D2[i, ]
            dim(D.mat) <- c(ncz, ncz)
            out[i] <- sum(D2[i, ] * cS.postVB, na.rm = TRUE) + sum((post.b %*% D.mat) * post.b, na.rm = TRUE)   
        }
        J <- jacobian2(attr(D, "L"), ncz)
        drop(0.5 * (n * D1 - out) %*% J)
    }
    c(score.y, score.t, score.b)
}
