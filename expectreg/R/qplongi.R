qplongi <-
function (design, lambda, types, p, y, P.list, Cmat, Bx, K, i.ja, 
    smooth, id) 
{
    n <- length(y)
    mp <- length(p)
    n_i <- table(id)
    n_obs <- length(unique(id))
    d <- length(design)
    Z_small <- matrix(0, ncol = n_obs, nrow = n)
    Z_small[1:n_i[1], 1] <- rep(1, n_i[1])
    for (j in 2:n_obs) {
        Z_small[(cumsum(n_i)[(j - 1)] + 1):(cumsum(n_i)[(j - 
            1)] + n_i[j]), j] <- rep(1, n_i[j])
    }
    Z_random <- rep(1, times = length(p)) %x% Z_small
    design[[d + 1]] <- list(B = Z_random, P = diag(n_obs), x = NULL, 
        type = "random", constraint = matrix(0, nrow = 2, ncol = ncol(Z_random)))
    design[[d + 1]]$xname <- "random intercept"
    types[[d + 1]] <- "random"
    K[d + 1] <- ncol(Z_random)
    P.list[[d + 1]] <- design[[d + 1]]$P
    Cmat[[length(design)]] = matrix(0, nrow = 2, ncol = n_obs)
    for (k in 2:length(p)) {
        Cmat[[length(design)]] = rbind(cbind(Cmat[[length(design)]], 
            matrix(0, nrow = nrow(Cmat[[length(design)]]), ncol = ncol(as.matrix(design[[length(design)]]$constraint)))), 
            cbind(matrix(0, nrow = nrow(as.matrix(design[[length(design)]]$constraint)), 
                ncol = ncol(Cmat[[length(design)]])), as.matrix(design[[length(design)]]$constraint)))
    }
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
        if (types[[i]] == "pspline" || types[[i]] == "random") {
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
        if (types != "random") 
            B <- Bp %x% Bx
        else {
            B <- Bx
        }
        W0 <- matrix(weights, nrow = nrow(B), ncol = ncol(B), 
            byrow = FALSE)
        W0B <- W0 * B
        if (mp > 1) {
            Amat <- B[-(1:n), ] - B[-((nrow(B) - n + 1):nrow(B)), 
                ]
            if (types == "intercept" || types == "random") {
                Amat <- matrix(0, ncol = ncol(Amat), nrow = nrow(Amat))
                meq <- 0
            }
            else meq <- mp
            if (!i.ja && i == 1) 
                meq <- 0
            if (ncol(Amat) == ncol(Cmat)) {
                Amat <- rbind(Amat, Cmat)
                bvec <- matrix(0, ncol = n * (mp - 1) + meq + 
                  nrow(Cmat), 1)
            }
            else bvec <- matrix(0, ncol = n * (mp - 1) + meq, 
                1)
            if ((types != "intercept" && i.ja && types != "random" && 
                types != "factor") || (types != "intercept" && 
                !i.ja && i != 1 && types != "random" && types != 
                "factor")) {
                Amat_pC_help <- matrix(0, ncol = n * mp, nrow = mp)
                for (m in 1:mp) Amat_pC_help[m, ((m - 1) * n + 
                  1):(m * n)] <- 1
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
                print("hier!")
            }
        }
        else {
            Amat <- matrix(0, nrow = ncol(Bx), ncol = ncol(Bx))
            bvec <- rep(0, ncol(Bx))
            meq <- 0
            W0B <- B
        }
        if (types != "random") 
            P_help <- diag(1, mp - 1 + bdegp, mp - 1 + bdegp) %x% 
                P
        else P_help <- P
        Dmat_part <- t(B) %*% W0B
        Dmat <- Dmat_part + lambda * P_help
        dvec <- t(yy) %*% W0B
        a_pq <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = t(Amat), 
            bvec = bvec, meq = meq)
        a_sol <- as.vector(unlist(a_pq$solution))
        z <- B %*% a_sol
        if (types == "random") 
            tau2dfi <- as.numeric((t(a_sol) %*% a_sol))
        else tau2dfi <- as.numeric((t(a_sol) %*% (P_help) %*% 
            a_sol))
        dfi <- sum(diag(solve(Dmat) %*% t(B) %*% W0B))
        ncexpect <- list(lambda = lambda, fitted = as.vector(z), 
            coefficients = a_sol, weights = W0[, 1], tau2dfi = tau2dfi, 
            Amat = Amat, bvec = bvec, dfi = dfi, Dmat_parti = Dmat_part, 
            P_helpi = P_help)
    }
    pw <- c()
    for (k in 1:mp) {
        pw <- c(pw, rep(p[k], n))
    }
    yy <- rep(y, mp)
    noncross.backfit <- function(yy, Bx, Bp, pw, bdegp, n, mp, 
        lambdanp, P.list, design, partbasis.list, Cmat, types, 
        weights = rep(0.5, times = n * mp)) {
        f <- matrix(0, nrow = n * mp, ncol = length(design))
        alpha <- mean(yy)
        ok <- TRUE
        rss0 <- 0
        coefficients <- list()
        coeff.vec <- list()
        Z <- list()
        part.resid <- list()
        Dmat_part <- list()
        P_help <- list()
        dfi <- vector(length = length(design))
        tau2dfi <- vector(length = length(design))
        count <- 0
        while (ok) {
            for (i in 1:length(design)) {
                difw <- 1
                it <- 0
                while ((difw > 0) && (it < 150)) {
                  weights0 <- weights
                  epsilon <- yy - rowSums(as.matrix(f[, -i])) - 
                    alpha
                  yywi <- rowSums(as.matrix(f[, -i])) + alpha
                  nc.obj <- noncross.fit(yy = epsilon, Bx = as.matrix(design[[i]][[1]]), 
                    Bp = Bp, mp = mp, bdegp = bdegp, lambda = lambdanp[i], 
                    P = as.matrix(P.list[[i]]), weights = weights, 
                    Cmat[[i]], types = types[[i]], yywh = yywi, 
                    i = i)
                  coeff.vec[[i]] <- nc.obj$coefficients
                  if (types[[i]] != "random") 
                    coefficients[[i]] <- matrix(nc.obj$coefficients, 
                      ncol = mp)
                  else coefficients[[i]] <- nc.obj$coefficients
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
                  ncol = mp)
                Z[[i]] <- matrix(f[, i], nrow = n, ncol = mp)
                dfi[i] <- nc.obj$dfi
                Dmat_part[[i]] <- nc.obj$Dmat_parti
                P_help[[i]] <- nc.obj$P_helpi
                if (types[[i]] != "intercept") 
                  dfi[i] <- dfi[i]
                tau2dfi[i] <- nc.obj$tau2dfi
            }
            rss <- sum((yy - yyhat)^2)
            if (abs(rss - rss0) < 1e-06 * rss) 
                ok <- FALSE
            rss0 <- rss
            count <- count + 1
            if (count == 150) 
                stop("No convergence within the backfitting algorithm. Stopping after 150 iterations.")
        }
        ncexpectbf <- list(lambdanp = lambdanp, fitted = f, coefficients = coefficients, 
            weights = weights, values = Z, alpha = alpha, yyhat = yyhat, 
            dfi = dfi, tau2dfi = tau2dfi, part.resid = part.resid, 
            coeff.vec = coeff.vec, Dmat_part = Dmat_part, P_help = P_help, 
            rss = rss)
    }
    tau2 <- rep(NA, length(lambdanp))
    sig2 <- NA
    if (smooth == "schall") {
        design_neu <- design
        design_neu[[(d + 1)]]$B <- Z_small
        for (i in 1:length(design)) {
            if ((types[[i]] == "parametric") || (types[[i]] == 
                "intercept")) 
                lambdanp[i] <- 0
            else if (types[[i]] == "pspline" || (types[[i]] == 
                "random")) {
                dc <- 1
                it <- 1
                lambda_einzel <- lambdanp[i]
                while (dc >= 0.001 && it < 100) {
                  print("hier")
                  mean.erg <- noncross.backfit(yy = y, Bx = Bx, 
                    Bp = c(1), pw = 0.5, bdegp = bdegp, n = n, 
                    mp = 1, lambdanp = lambdanp, P.list = P.list, 
                    design = design_neu, partbasis.list = partbasis.list, 
                    Cmat = Cmat, types = types)
                  tau2[i] <- (mean.erg$tau2dfi[i]/mean.erg$dfi[i]) + 
                    1e-06
                  print(paste("i =", i))
                  print("tau2[i]")
                  print(tau2[i])
                  print("mean.erg$dfi")
                  print(mean.erg$dfi[i])
                  sig2df <- mean.erg$rss
                  sig2 <- as.numeric(sig2df/n)
                  print(paste("n =", n))
                  print("sig2")
                  print(sig2)
                  lambda_neu <- sig2/tau2[i]
                  dc <- abs(log10(lambda_einzel) - log10(lambda_neu))
                  print("dc")
                  print(dc)
                  print("lambda_alt")
                  print(lambda_einzel)
                  print("lambda_neu")
                  print(lambda_neu)
                  it <- it + 1
                  lambda_einzel <- lambda_neu
                  lambdanp[i] <- lambda_neu
                  if (it == 100) 
                    warning("Schall algorithm did not converge. Stopping after 100 iterations.")
                }
                lambdanp[i] <- lambda_einzel
            }
        }
        nc.bf <- noncross.backfit(yy = yy, Bx = Bx, Bp = Bp, 
            pw = pw, bdegp = bdegp, n = n, mp = mp, lambdanp = lambdanp, 
            P.list = P.list, design = design, partbasis.list = partbasis.list, 
            Cmat = Cmat, types = types)
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
            Cmat = Cmat, types = types)
    }
    else {
        nc.bf <- noncross.backfit(yy = yy, Bx = Bx, Bp = Bp, 
            pw = pw, bdegp = bdegp, n = n, mp = mp, lambdanp = lambdanp, 
            P.list = P.list, design = design, partbasis.list = partbasis.list, 
            Cmat = Cmat, types = types)
    }
    nc.bf$tau2 = tau2
    nc.bf$sig2 = sig2
    nc.bf
}
