Score.splineGH <-
function (thetas) {
    thetas <- relist(thetas, skeleton = list.thetas)
    betas <- thetas$betas
    sigma <- exp(thetas$log.sigma)
    gammas <- thetas$gammas
    gammas.bs <- thetas$gammas.bs
    alpha <- thetas$alpha
    Dalpha <- thetas$Dalpha
    D <- thetas$D
    D <- if (diag.D) exp(D) else chol.transf(D)
    eta.yx <- as.vector(X %*% betas)
    eta.tw1 <- if (!is.null(W1)) as.vector(W1 %*% gammas) else rep(0, n)
    eta.tw2 <- as.vector(W2 %*% gammas.bs)
    if (parameterization %in% c("value", "both")) {
        Y <- as.vector(Xtime %*% betas) + Ztime.b
        Ys <- as.vector(Xs %*% betas) + Zsb
        WintF.vl.alph <- c(WintF.vl %*% alpha)
        Ws.intF.vl.alph <- c(Ws.intF.vl %*% alpha)
        eta.t <- eta.tw2 + eta.tw1 + WintF.vl.alph * Y
        eta.s <- Ws.intF.vl.alph * Ys
    }
    if (parameterization %in% c("slope", "both")) {
        Y.deriv <- as.vector(Xtime.deriv %*% betas[indFixed]) + Ztime.b.deriv
        Ys.deriv <- as.vector(Xs.deriv %*% betas[indFixed]) + Zsb.deriv
        WintF.sl.alph <- c(WintF.sl %*% Dalpha)
        Ws.intF.sl.alph <- c(Ws.intF.sl %*% Dalpha)
        eta.t <- if (parameterization == "both")
            eta.t + WintF.sl.alph * Y.deriv
        else
            eta.tw2 + eta.tw1 + WintF.sl.alph * Y.deriv
        eta.s <- if (parameterization == "both")
            eta.s + Ws.intF.sl.alph * Ys.deriv
        else
            Ws.intF.sl.alph * Ys.deriv
    }    
    eta.ws <- as.vector(W2s %*% gammas.bs)
    exp.eta.tw.P <- exp(eta.tw1) * P
    mu.y <- eta.yx + Ztb
    logNorm <- dnorm(y, mu.y, sigma, TRUE)
    log.p.yb <- rowsum(logNorm, id)
    log.hazard <- eta.t
    Int <- wk * exp(eta.ws + eta.s)
    log.survival <- - exp.eta.tw.P * rowsum(Int, id.GK, reorder = FALSE)
    dimnames(log.survival) <- NULL
    log.p.tb <- rowsum(d * log.hazard + log.survival, idT, reorder = FALSE)
    log.p.b <- if (control$typeGH == "simple") {
        rep(dmvnorm(b, rep(0, ncz), D, TRUE), each = n)
    } else {
        matrix(dmvnorm(do.call(rbind, lis.b), rep(0, ncz), D, TRUE), 
            n, k, byrow = TRUE)
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
    sc2 <- numeric(ncx)
    for (i in 1:ncx) {
        ki <- exp.eta.tw.P * switch(parameterization,
            "value" = rowsum(Int * Ws.intF.vl.alph * Xs[, i], id.GK, reorder = FALSE),
            "slope" = {ii <- match(i, indFixed); 
                if (is.na(ii)) 0 else rowsum(Int * Ws.intF.sl.alph * Xs.deriv[, ii], 
                    id.GK, reorder = FALSE)},
            "both" = {ii <- match(i, indFixed);
                rowsum(Int * (Ws.intF.vl.alph * Xs[, i] + 
                    Ws.intF.sl.alph * if (is.na(ii)) 0 else Xs.deriv[, ii]), 
                        id.GK, reorder = FALSE)}
        )
        ki <- c(rowsum(ki, idT, reorder = FALSE))
        kii <- c((p.byt * ki) %*% wGH)
        sc2[i] <- switch(parameterization,
            "value" = {
                ddd <- tapply(d * WintF.vl.alph * Xtime[, i], idT, sum)
                - sum(ddd - kii, na.rm = TRUE)
            },
            "slope" = {
                ii <- match(i, indFixed)
                if (is.na(ii)) 0 else {
                    ddd <- tapply(d * WintF.sl.alph * Xtime.deriv[, ii], idT, sum)
                    - sum(ddd - kii, na.rm = TRUE)
                }
            },
            "both" = {
                ii <- match(i, indFixed)
                ddd <- tapply(d * (WintF.vl.alph * Xtime[, i] + WintF.sl.alph * 
                    if (is.na(ii)) 0 else Xtime.deriv[, ii]), idT, sum)
                - sum(ddd - kii, na.rm = TRUE)
            })
    }    
    score.y <- c(sc1 + sc2, - sigma * (- N / sigma + 
        drop(crossprod(mu, mu - 2 * Zb) + crossprod(Zb) + tr.tZZvarb) / sigma^3))
    scgammas1 <- if (!is.null(W1)) {
        scg1 <- numeric(ncol(W1))
        for (jj in seq_along(scg1)) {
            tt <- rowsum(W1[, jj] * log.survival, idT, reorder = FALSE)
            scg1[jj] <- sum(c((p.byt * tt) %*% wGH), na.rm = TRUE)
        }
        - colSums(W1 * d, na.rm = TRUE) - scg1
    } else 
        NULL
    scgammas2 <- numeric(nk)
    for (i in 1:nk) {
        kk <- exp.eta.tw.P * rowsum(Int * W2s[, i], id.GK, reorder = FALSE)
        kk <- rowsum(kk, idT, reorder = FALSE)
        scgammas2[i] <- - sum(W2[, i] * d) + sum(c((p.byt * kk) %*% wGH))
    }
    scalpha <- if (parameterization %in% c("value", "both")) {
        rr <- numeric(ncol(WintF.vl))
        for (l in seq_along(rr)) {
            rrr <- exp.eta.tw.P * rowsum(Int * Ws.intF.vl[, l] * 
                    Ys, id.GK, reorder = FALSE)
            rrr <- rowsum(rrr, idT, reorder = FALSE)
            rr[l] <- - sum((p.byt * (rowsum(d * WintF.vl[, l] * Y, idT, 
                    reorder = FALSE) - rrr)) %*% wGH, na.rm = TRUE)
        }
        rr
    } else NULL
    scalpha.D <- if (parameterization %in% c("slope", "both")) {
        rr <- numeric(ncol(WintF.sl))
        for (l in seq_along(rr)) {
            rrr <- exp.eta.tw.P * rowsum(Int * Ws.intF.sl[, l] * 
                    Ys.deriv, id.GK, reorder = FALSE)
            rrr <- rowsum(rrr, idT, reorder = FALSE)            
            rr[l] <- - sum((p.byt * (rowsum(d * WintF.sl[, l] * Y.deriv, 
                    idT, reorder = FALSE) - rrr)) %*% wGH, na.rm = TRUE)
        }
        rr
    } else NULL
    score.t <- c(scgammas1, scalpha, scalpha.D, scgammas2)
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
