fitModel <- function (k1, k2, Y, n.surr){
    Y <- as.matrix(Y)
    G <- nrow(Y)
    n <- k1 + k2
    if (n.surr == 0) {
        SSE <- 0
        for (i in 1:G) {
            for (u in 1:n) {
                SSE <- SSE + ifelse(u <= k1, (Y[i, u] - mean(Y[i, 
                  1:k1]))^2, (Y[i, u] - mean(Y[i, (k1 + 1):n]))^2)
            }
        }
        MLE <- SSE/(n*G)
        MSE <- SSE/((n - 2) * G)
        G.hat <- V1.hat <- GV1.hat <- rep(0, G)
        for (i in 1:G) {
            G.hat[i] <- mean(Y[i, ]) - mean(Y)
            V1.hat[i] <- mean(Y[, 1:k1]) - mean(Y)
            GV1.hat[i] <- mean(Y[i,1:k1]) - mean(Y[i,]) - mean(Y[,1:k1]) + mean(Y)
        }
        mu.hat <- mean(Y)
        V2.hat <- -(k1/k2) * V1.hat
        GV2.hat <- -(k1/k2) * GV1.hat
        V.hat <- c(V1.hat, V2.hat)
        GV.hat <- cbind(GV1.hat, GV2.hat)
        AIC <- n * G * (log(2 * pi) + log(MLE)) + n*G + 2*(2*G + 1)
        coef1 <- list(mu.hat = mu.hat, G.hat = G.hat, V.hat = V.hat, 
            GV.hat = GV.hat, MSE = MSE, AIC = AIC)
    }
    if (n.surr > 0) {
        m_lab1 <- matrix(rep(apply(Y[, 1:k1], 1, mean), k1), 
            G, k1, byrow = F)
        m_lab2 <- matrix(rep(apply(Y[, (k1 + 1):n], 1, mean), 
            k2), G, k2, byrow = F)
        e <- Y - cbind(m_lab1, m_lab2)
        pls <- mvr(t(e) ~ t(Y), ncomp = n.surr, method = "oscorespls")
        sc <- scores(pls)
        U <- var(sc[1:k1, 1]) * (k1 - 1)/k1
        Q <- matrix(0, n.surr, n.surr)
        B2 <- V <- E <- H <- N <- N1 <- K <- P <- h <- rep(0, 
            n.surr)
        Z1.bar <- Z2.bar <- S1 <- U1 <- C2 <- S2 <- L <- rep(0, 
            n.surr)
        Z1 <- Z2 <- Z <- vec <- T <- set <- list()
        T1 <- T2 <- matrix(0, 2, n)
        C.mat <- rep(0, n * G)
        for (u in 1:n.surr) {
            N[u] <- (k1/n) * (mean(sc[1:k1, u]) - mean(sc[(k1 + 
                1):n, u]))
            N1[u] <- (k1/n) * (mean(sc[1:k1, u]^2) - mean(sc[(k1 + 
                1):n, u]^2)) + N[u]
            K[u] <- (N[u] * (N[u] + mean(sc[1:k1, u])) - N1[u])/mean(sc[, 
                u]^2)
            h[u] <- mean(sc[, u]) - mean(sc[1:k1, u])
            P[u] <- k1/n * (cov(sc[(k1 + 1):n, 1], sc[(k1 + 1):n, 
                u]) * (k2 - 1)/k2 - cov(sc[1:k1, 1], sc[1:k1, 
                u]) * (k1 - 1)/k1)
            Z[[u]] <- matrix(0, G, n)
            Z1[[u]] <- matrix(0, G, k1)
            Z2[[u]] <- matrix(0, G, k2)
            Z[[u]] <- matrix(rep(sc[, u], G), G, n, byrow = TRUE)
            Z1[[u]] <- matrix(rep(sc[1:k1, u], G), G, k1, byrow = TRUE)
            Z2[[u]] <- matrix(rep(sc[(k1 + 1):n, u], G), G, k2, 
                byrow = TRUE)
            Z1.bar[u] <- mean(sc[1:k1, u])
            Z2.bar[u] <- mean(sc[(k1 + 1):n, u])
            B2[u] <- (cov(sc[1:k1, u], sc[1:k1, 1]) * (k1 - 1)/k1)/((k1 - 
                1)/k1 * var(sc[1:k1, 1]))
            S2[u] <- var(sc[, u]) * (n - 1)/n
            L[u] <- (k1/n) * (mean(sc[1:k1, u]^2) - mean(sc[(k1 + 
                1):n, u]^2))
            U1[u] <- N[u] * (mean(sc[, u]) + mean(sc[1:k1, u])) - 
                N[u]^2 - L[u]
            C2[u] <- cov(sc[, u], sc[, 1]) * (n - 1)/n + U1[1] * 
                B2[u] + N[1] * h[u]
            F[u] <- S2[u] + N[u] * h[u]
            A <- rep(sc[, u]/(n * G) - mean(sc[, u])/(n * G) + 
                N[u]/(n * G), G)
            ce <- (N[u] * (N[1] - mean(sc[1:k1, 1])) - P[u])/(U * 
                k1 * G)
            B <- -N[u]/(G * k1) - ce * (sc[1:k1, 1] - mean(sc[1:k1, 
                1]))
            B <- rep(c(B, rep(0, k2)), G)
            C.mat <- rbind(C.mat, A + B)
        }
        C.mat <- C.mat[-1, ]
        for (u in 1:n.surr) {
            for (l in 1:n.surr) {
                Q[l, u] <- cov(sc[, u], sc[, l]) * (n - 1)/n + 
                  (P[l] - (N[1] - mean(sc[1:k1, 1])) * N[l]) * 
                    B2[u] + N[l] * h[u]
            }
        }
        set[[1]] <- 1:k1
        set[[2]] <- (k1 + 1):n
        for (j in 1:2) {
            for (k in set[[j]]) {
                T1[j, k] <- (sc[k, 1] - mean(sc[, 1]))
                T2[j, k] <- U1[1]/(U * k1 * G) * (sc[k, 1] - 
                  mean(sc[1:k1, 1])) - N[1]/(G * k1)
            }
        }
        T1 <- c(T1[1, set[[1]]], T1[2, set[[2]]])
        sc.means <- as.vector(apply(sc, 2, mean))
        B1 <- (cov(as.vector(Y[, 1:k1]), as.vector(Z1[[1]])) * 
            (G * k1 - 1)/(G * k1))/(var(sc[1:k1, 1]) * (k1 - 
            1)/k1)
        T.vals <- solve(Q) %*% C.mat
        for (u in 1:n.surr) T[[u]] <- matrix(T.vals[u, ], G, 
            n, byrow = TRUE)
        beta.hat <- as.vector(T.vals %*% as.matrix(as.vector(t(Y))))
        VZ1.hat <- B1 - sum(B2 * beta.hat)
        VZ2.hat <- (-k1/k2) * VZ1.hat
        mu.hat <- mean(Y) - sum(beta.hat * sc.means) - N[1] * 
            VZ1.hat
        V1.hat <- mean(Y[, 1:k1]) - mu.hat - sum(beta.hat * Z1.bar) - 
            VZ1.hat * Z1.bar[1]
        V2.hat <- -(k1/k2) * V1.hat
        V.A <- (S2[1]/F[1])^2 * (1 - 1/G) * (1/k1 - 1/n)
        p1 <- sum(T1^2/n^2 + 2 * N[1]/(n^2 * G) * T1) + N[1]^2/(n * 
            G)
        p2 <- sum(G * T2[1, 1:k1]^2 + (2/n) * T1[1:k1] * T2[1, 
            1:k1] + (2 * N[1]/n) * T2[1, 1:k1])
        V.B <- (h[1]/F[1])^2 * (p1 + p2)
        cnt <- 0
        for (u1 in 1:n.surr) {
            for (u2 in 1:n.surr) {
                cnt <- cnt + C2[u1] * C2[u2] * sum(T[[u1]] * 
                  T[[u2]])
            }
        }
        V.C <- (h[1]/F[1])^2 * cnt
        C.AB <- (1/n) * (1 - 1/G) * (S2[1] * h[1]/F[1]^2) * (mean(T1[1:k1]) - 
            mean(T1))
        G.hat <- GV1.hat <- GV2.hat <- GZ1.hat <- val <- rep(0, 
            G)
        vhat.gv1h <- vhat.gv2h <- rep(0, G)
        for (i in 1:G) {
            val[i] <- mean(Y[i, 1:k1]) - mean(Y[i, ]) - mean(Y[, 
                1:k1]) + mean(Y)
            P1 <- S2[1]/F[1] * val[i]
            P2 <- h[1]/F[1] * (cov(Y[i, ], sc[, 1]) * (n - 1)/n + 
                U1[1] * B1 - N[1] * (mean(Y[, 1:k1]) - mean(Y)))
            cnt1 <- cnt2 <- cnt3 <- 0
            for (u in 1:n.surr) {
                cnt1 <- cnt1 + C2[u] * (mean(T[[u]][i, 1:k1]) - 
                  mean(T[[u]][i, ]) - mean(T[[u]][, 1:k1]) + 
                  mean(T[[u]]))
                cnt2 <- cnt2 + C2[u] * (mean(T1 * T[[u]][i, ]) + 
                  sum(T2[1, 1:k1] * t(T[[u]][, 1:k1])) + N[1] * 
                  mean(T[[u]]))
                cnt3 <- cnt3 + C2[u] * sum(T[[u]] * Y)
            }
            P3 <- -(h[1]/F[1]) * cnt3
            GV1.hat[i] <- P1 + P2 + P3
            GV2.hat[i] <- (-k1/k2) * GV1.hat[i]
            GZ1.hat[i] <- (GV1.hat[i] - val[i])/h[1]
            G.hat[i] <- mean(Y[i, ]) - mu.hat - sum(beta.hat * 
                apply(sc, 2, mean)) - GZ1.hat[i] * mean(sc[, 
                1]) - N[1] * VZ1.hat
            C.AC <- -(S2[1] * h[1]/F[1]^2) * cnt1
            C.BC <- -(h[1]/F[1])^2 * cnt2
            vhat.gv1h[i] <- V.A + V.B + V.C + 2 * C.AB + 2 * 
                C.AC + 2 * C.BC
            vhat.gv2h[i] <- (k1/k2)^2 * vhat.gv1h[i]
        }
        V.hat <- c(V1.hat, V2.hat)
        GV.hat <- cbind(GV1.hat, GV2.hat)
        VZ1.hat <- c(VZ1.hat, VZ2.hat)
        vhat.gvh <- cbind(vhat.gv1h, vhat.gv2h)
        A1 <- Y[, 1:k1] - matrix(rep(apply(Y[, 1:k1], 1, mean), 
            each = k1), G, k1, byrow = TRUE)
        A2 <- Y[, (k1 + 1):n] - matrix(rep(apply(Y[, (k1 + 1):n], 
            1, mean), each = k2), G, k2, byrow = TRUE)
        W1 <- matrix(rep(apply(beta.hat * (Z1.bar - t(sc[1:k1, 
            ])), 2, sum), G), G, k1, byrow = TRUE)
        W2 <- matrix(rep(apply(beta.hat * (Z2.bar - t(sc[(k1 + 
            1):n, ])), 2, sum), G), G, k2, byrow = TRUE)
        P1 <- GZ1.hat * (mean(sc[1:k1, 1]) - Z1[[1]])
        P2 <- GZ1.hat * (mean(sc[(k1 + 1):n, 1]) - Z2[[1]])
        D1 <- matrix(rep(VZ1.hat[1] * (mean(sc[1:k1, 1]) - sc[1:k1, 
            1]), G), G, k1, byrow = TRUE)
        D2 <- matrix(rep(VZ1.hat[2] * (mean(sc[(k1 + 1):n, 1]) - 
            sc[(k1 + 1):n, 1]), G), G, k2, byrow = TRUE)
        A <- cbind(A1, A2)
        W <- cbind(W1, W2)
        P <- cbind(P1, P2)
        D <- cbind(D1, D2)
        ERR <- A + W + P + D
        SSE <- sum(ERR^2)
        MLE <- SSE/(n*G)
        MSE <- SSE/(n * G - 3 * G - n.surr)
        AIC <- n * G * (log(2 * pi) + log(MLE)) + n*G + 2*(3*G + n.surr + 1) 
        coef2 <- list(mu.hat = mu.hat, G.hat = G.hat, V.hat = V.hat, 
            GV.hat = GV.hat, sc = sc, beta.hat = beta.hat, GZ1.hat = GZ1.hat, 
            VZ1.hat = VZ1.hat, vhat.gvh = vhat.gvh, MSE = MSE, 
            AIC = AIC)
    }
    if (n.surr == 0) 
        coef <- coef1
    if (n.surr != 0) 
        coef <- coef2
    class(coef) <- c("fitModel", "list", "vector")
    return(coef)
}

## new summary function S3
summary.fitModel <- function(object){
cat("mu.hat: \n", object$mu.hat, "\n")
cat("G.hat: \n")
print(object$G.hat) 
cat("\n")
cat("V.hat: \n")
print(object$V.hat)
cat("\n")
cat("GV.hat: \n")
print(object$GV.hat)
cat("\n")
cat("Z: \n")
print(object$sc)
cat("\n")
cat("beta.hat: \n")
print(object$beta.hat)
cat("\n")
cat("GZ1.hat: \n")
print(object$GZ1.hat)
cat("\n")
cat("VZ1.hat: \n")
print(object$VZ1.hat)
cat("\n")
cat("Vhat.GV.hat: \n")
print(object$vhat.gvh)
cat("\n")
cat("MSE: \n")
print(object$MSE)
cat("\n")
cat("AIC: \n")
print(object$AIC)
cat("\n")
}

## new print function S3
print.fitModel <- function(x){
cat("Estimated coefficients of the surrogate variables: \n")
print(x$beta.hat)
cat("\n")
cat("Estimated Mean Squared Error of the fitted model: \n")
print(x$MSE)
cat("\n")
cat("AIC value of the fitted model: \n")
print(x$AIC)
cat("\n")
}