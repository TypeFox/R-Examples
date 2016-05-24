qpvzentr <-
function (design, lambda, types, p, y, P.list, Cmat, Bx, K, i.ja, 
    smooth) 
{
    n <- length(y)
    mp <- length(p)
    if (length(lambda) < length(design)) {
        lambda <- rep(lambda[1], length(design))
    }
    bdegp <- 1
    p.knots = p
    p.knots[1] <- p.knots[1] - 1e-07
    p.knots[mp] <- p.knots[mp] + 1e-07
    dx.left <- ((p.knots[1] + 1/(mp + 1))/bdegp)
    outer.left <- seq(p.knots[1] - bdegp * dx.left, p.knots[1] - 
        dx.left, by = dx.left)
    dx.right <- (((1 + 1/(mp + 1)) - p.knots[mp])/bdegp)
    outer.right <- seq(p.knots[mp] + dx.right, p.knots[mp] + 
        bdegp * dx.right, by = dx.right)
    knots <- c(outer.left, p.knots, outer.right)
    Bp <- splineDesign(knots, p.knots, bdegp + 1)
    dfparam <- vector(length = length(types))
    for (i in 1:length(types)) {
        dfparam[i] <- (types[[i]] == "parametric")
    }
    dfparam <- sum(dfparam) + i.ja
    bigBxnp <- matrix(0, ncol = 0, nrow = n)
    bigPnp <- matrix(0, ncol = 0, nrow = (21 - 2))
    lambdanp <- rep(0, times = length(lambda))
    indexnp <- vector(length = 0)
    for (i in 1:length(types)) {
        if (types[[i]] == "pspline") {
            bigBxnp <- cbind(bigBxnp, design[[i]][[1]])
            bigPnp <- rbind(cbind(bigPnp, matrix(0, nrow = nrow(bigPnp), 
                ncol = ncol(design[[i]][[2]]))), cbind(matrix(0, 
                nrow = nrow(design[[i]][[2]]), ncol = ncol(bigPnp)), 
                design[[i]][[2]]))
            lambdanp[i] <- lambda[i]
            indexnp <- c(indexnp, i)
        }
    }
    partbasis.list = list()
    for (k in 1:length(design)) {
        partbasis.list[[k]] = (sum(K[0:(k - 1)]) + 1):(sum(K[0:k]))
    }
    acv.noncross <- function(penalty, yy, Bx, Bp, pw, bdegp, 
        n, mp, P.list, design, partbasis.list) {
        aa <- noncross.backfit(yy, Bx, Bp, pw, bdegp, n, mp, 
            abs(penalty), P.list, design, partbasis.list, Cmat = Cmat)
        score = (aa$weights * (yy - aa$yyhat)^2)/((1 - sum(aa$dfi)/(n * 
            mp))^2)
        mean(score[which(is.finite(score))], na.rm = TRUE)
    }
    noncross.fit <- function(yy, Bx, Bp, bdegp, mp, lambda, P, 
        weights, Cmat, types, yywh, i) {
        B <- Bp %x% Bx
        W0 <- matrix(weights, nrow = nrow(B), ncol = ncol(B), 
            byrow = FALSE)
        W0B <- W0 * B
        Amat <- B[-(1:n), ] - B[-((nrow(B) - n + 1):nrow(B)), 
            ]
        if (types == "intercept") {
            Amat <- matrix(0, ncol = ncol(Amat), nrow = nrow(Amat))
            meq <- 0
        }
        else meq <- mp
        if (types != "intercept" && !i.ja && i == 1) 
            meq <- 0
        if (ncol(Amat) == ncol(Cmat)) {
            Amat <- rbind(Amat, Cmat)
            bvec <- matrix(0, ncol = n * (mp - 1) + meq + nrow(Cmat), 
                1)
        }
        else bvec <- matrix(0, ncol = n * (mp - 1) + meq, 1)
        if ((types != "intercept" && i.ja) || (types != "intercept" && 
            !i.ja && i != 1)) {
            Amat_pC_help <- matrix(0, ncol = n * mp, nrow = mp)
            for (m in 1:mp) Amat_pC_help[m, ((m - 1) * n + 1):(m * 
                n)] <- 1
            Amat_partCenter <- Amat_pC_help %*% B
            Amat <- rbind(Amat_partCenter, Amat)
            bvec_partNC <- yywh[-((nrow(B) - n + 1):nrow(B))] - 
                yywh[-(1:n)]
            bvec[, (meq + 1):(n * (mp - 1) + meq)] <- bvec_partNC
        }
        if (types != "intercept" && !i.ja && i == 1) {
            bvec_partNC <- yywh[-((nrow(B) - n + 1):nrow(B))] - 
                yywh[-(1:n)]
            bvec[, (meq + 1):(n * (mp - 1) + meq)] <- bvec_partNC
        }
        Dmat <- t(B) %*% W0B + 2 * lambda * (diag(1, mp - 1 + 
            bdegp, mp - 1 + bdegp) %x% P)
        dvec <- t(yy) %*% W0B
        a_pq <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = t(Amat), 
            bvec = bvec, meq = meq)
        a_sol <- as.vector(unlist(a_pq$solution))
        z <- B %*% a_sol
        tau2dfi <- as.numeric((t(a_sol) %*% (diag(1, mp - 1 + 
            bdegp, mp - 1 + bdegp) %x% P) %*% a_sol))
        dfi <- sum(diag(solve(Dmat) %*% t(B) %*% W0B))
        ncexpect <- list(lambda = lambda, fitted = as.vector(z), 
            coefficients = a_sol, weights = W0[, 1], tau2dfi = tau2dfi, 
            Amat = Amat, bvec = bvec, dfi = dfi)
    }
    pw <- c()
    for (k in 1:mp) {
        pw <- c(pw, rep(p[k], n))
    }
    yy <- rep(y, mp)
    noncross.backfit <- function(yy, Bx, Bp, pw, bdegp, n, mp, 
        lambdanp, P.list, design, partbasis.list, Cmat) {
        weights <- rep(0.5, times = n * mp)
        f <- matrix(0, nrow = n * mp, ncol = length(design))
        alpha <- mean(yy)
        ok <- TRUE
        rss0 <- 0
        coefficients <- list()
        coeff.vec <- list()
        Z <- list()
        part.resid <- list()
        dfi <- vector(length = length(design))
        tau2dfi <- vector(length = length(design))
        while (ok) {
            for (i in 1:length(design)) {
                difw <- 1
                it <- 0
                while ((difw > 0) && (it < 150)) {
                  weights0 <- weights
                  epsilon <- yy - rowSums(as.matrix(f[, -i])) - 
                    alpha
                  yywi <- rowSums(as.matrix(f[, -i])) + alpha
                  nc.obj <- noncross.fit(yy = epsilon, Bx = as.matrix(Bx[, 
                    partbasis.list[[i]]]), Bp = Bp, mp = mp, 
                    bdegp = bdegp, lambda = lambdanp[i], P = as.matrix(P.list[[i]]), 
                    weights = weights, Cmat[[i]], types = types[[i]], 
                    yywh = yywi, i = i)
                  coeff.vec[[i]] <- nc.obj$coefficients
                  coefficients[[i]] <- matrix(nc.obj$coefficients, 
                    ncol = length(p))
                  f[, i] <- nc.obj$fitted
                  yyhat <- alpha + rowSums(f)
                  weights <- pw * (yy > yyhat) + (1 - pw) * (yy <= 
                    yyhat)
                  it <- it + 1
                  difw <- sum(abs(weights0 - weights))
                  if (it == 150) 
                    warning("Weights did not converge. Stopping after 150 iterations.")
                }
                part.resid[[i]] <- matrix(epsilon, nrow = n, 
                  ncol = length(p))
                Z[[i]] <- matrix(f[, i], nrow = n, ncol = length(p))
                dfi[i] <- nc.obj$dfi
                if (types[[i]] != "intercept") 
                  dfi[i] <- dfi[i] - 1
                tau2dfi[i] <- nc.obj$tau2dfi
            }
            rss <- sum((yy - yyhat)^2)
            if (abs(rss - rss0) < 1e-06 * rss) 
                ok <- FALSE
            rss0 <- rss
        }
        ncexpectbf <- list(lambdanp = lambdanp, fitted = f, coefficients = coefficients, 
            weights = weights, values = Z, alpha = alpha, yyhat = yyhat, 
            dfi = dfi, tau2dfi = tau2dfi, part.resid = part.resid, 
            coeff.vec = coeff.vec)
    }
    if (smooth == "schall") {
        dc <- 1
        it <- 1
        while (dc >= 1e-04 && it < 100) {
            nc.bf <- noncross.backfit(yy = yy, Bx = Bx, Bp = Bp, 
                pw = pw, bdegp = bdegp, n = n, mp = mp, lambdanp = lambdanp, 
                P.list = P.list, design = design, partbasis.list = partbasis.list, 
                Cmat = Cmat)
            tau2 <- (nc.bf$tau2dfi/nc.bf$dfi) + 1e-06
            sig2df <- as.numeric(t(yy - nc.bf$yyhat) %*% (nc.bf$weights * 
                (yy - nc.bf$yyhat)))
            df <- sum(nc.bf$dfi)
            sig2 <- as.numeric(sig2df/(length(yy) - df))
            lambdanp <- sig2/tau2
            for (i in 1:length(design)) {
                if ((types[[i]] == "parametric") || (types[[i]] == 
                  "intercept")) 
                  lambdanp[i] <- 0
            }
            dc <- sum((nc.bf$lambdanp - lambdanp)^2)
            it <- it + 1
            if (it == 100) 
                warning("Schall algorithm did not converge. Stopping after 100 iterations.")
        }
    }
    else if (smooth == "acv") {
        acv.min <- nlm(acv.noncross, p = lambdanp, yy = yy, Bx = Bx, 
            Bp = Bp, pw = pw, bdegp = bdegp, n = n, mp = mp, 
            P = P.list, design = design, partbasis.list = partbasis.list, 
            ndigit = 8, iterlim = 50, gradtol = 1e-04)
        lambdanp <- abs(acv.min$estimate)
        nc.bf <- noncross.backfit(yy = yy, Bx = Bx, Bp = Bp, 
            pw = pw, bdegp = bdegp, n = n, mp = mp, lambdanp = lambdanp, 
            P.list = P.list, design = design, partbasis.list = partbasis.list, 
            Cmat = Cmat)
    }
    else {
        nc.bf <- noncross.backfit(yy = yy, Bx = Bx, Bp = Bp, 
            pw = pw, bdegp = bdegp, n = n, mp = mp, lambdanp = lambdanp, 
            P.list = P.list, design = design, partbasis.list = partbasis.list, 
            Cmat = Cmat)
        tau2 = NA
        sig2 = NA
    }
    nc.bf$tau2 = tau2
    nc.bf$sig2 = sig2
    nc.bf
}
