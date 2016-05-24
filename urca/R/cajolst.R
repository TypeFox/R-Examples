##
## cajolst
##
cajolst <- function (x, trend = TRUE, K = 2, season = NULL) 
{
    x <- as.matrix(x)
    K <- as.integer(K)
    if(K < 2){
      stop("\nK must be at least K=2.\n")
    }
    P <- ncol(x)
    arrsel <- P
    N <- nrow(x)
    if (!is.null(season)) {
        s <- season - 1
    }
    else {
        s <- 0
    }
    if (N * P < P + s * P + K * P^2 + P * (P + 1)/2) 
        stop("\nInsufficient degrees of freedom.\n")
    if (P > 5) 
        warning("\nToo many variables, critical values cannot be computed.\n")
    if (!(is.null(season))) {
        dum <- (diag(season) - 1/season)[, -season]
        dums <- dum
        while (nrow(dums) < N) {
            dums <- rbind(dums, dum)
        }
        dums <- dums[1:N, ]
        if (NA %in% x) {
            idx.NA <- 1:N
            ind <- as.logical(sapply(idx.NA, function(z) sum(is.na(x[z, 
                ]) * 1)))
            ind2 <- ind * (1:N)
            dums <- dums[-ind2, ]
        }
    }
    x2 <- na.omit(x)
    Ntot <- nrow(x2)
    y <- embed(x2, (K + 1))
    rhs <- y[, -c(1:P)]
    if (!trend) {
        rhs <- y[, -c(1:P)]
    }
    else {
        trd <- seq(K + 1, nrow(y) + K)
        rhs <- cbind(trd, y[, -c(1:P)])
    }
    N <- nrow(y)
    if (!(is.null(season))) {
        rhs <- cbind(dums[-(1:K), ], rhs)
    }
    lhs <- y[, 1:P]
    idx <- 1:(N - 1)
    tau <- function(t) {
        dt <- c(rep(0, t), rep(1, N - t))
        det(crossprod(resid(lm(lhs ~ dt + rhs))))
    }
    tau.hat <- sapply(idx, tau)
    tau.opt <- which.min(tau.hat) + K
    tau.bp <- tau.opt + 1
    dt <- c(rep(0, tau.opt), rep(1, N - tau.opt))
    if(!trend & is.null(season)){
      rhs.aux <- dt
    } else {
      rhs.aux <- cbind(dt, rhs[, -c((ncol(rhs)-K*ncol(x)+1):ncol(rhs))])
    }
    reg.opt <- lm(lhs ~ rhs.aux)
    dt <- c(rep(0, tau.opt), rep(1, Ntot - tau.opt))
    uv <- c(rep(1, Ntot))
    if (!trend) {
        if (!is.null(season)) {
            yfit <- x - uv%*%t(coef(reg.opt)[1, ]) - dt %*% t(coef(reg.opt)[2, ]) - dums %*% coef(reg.opt)[3:(2 + season - 1), ]
        }else{
            yfit <- x - uv%*%t(coef(reg.opt)[1, ]) - dt %*% t(coef(reg.opt)[2, ])
        }
    }else if (trend){
        trd <- 1:Ntot
        if (!is.null(season)) {
            yfit <- x - uv%*%t(coef(reg.opt)[1, ]) - dt %*% t(coef(reg.opt)[2, ]) - dums %*% coef(reg.opt)[3:(2 + season - 1), ] - trd %*% t(coef(reg.opt)[season + 2, ])
        }else{
            yfit <- x - uv%*%t(coef(reg.opt)[1, ]) - dt %*% t(coef(reg.opt)[2, ]) - trd %*% t(coef(reg.opt)[3, ])
        }
    }
    x <- na.omit(yfit)
    N <- nrow(x)
    spec <- "transitory"
    Z <- embed(diff(x), K)
    Z0 <- Z[, 1:P]
    Z1 <- Z[, -c(1:P)]
    ZK <- x[-N, ][K:(N - 1), ]
    idx <- 0:(P - 1)
    if (trend) {
      cvals <- matrix(c(5.423, 13.784, 25.931, 42.083, 61.918, 6.785, 15.826, 28.455, 45.204, 65.662, 10.042, 19.854, 33.757, 51.601, 73.116), nrow=5, ncol=3)
      model <- "with linear trend in shift correction"
    }else if(!trend){
      cvals <- matrix(c(2.996, 10.446, 21.801, 36.903, 55.952, 4.118, 12.276, 24.282, 40.067, 59.749, 6.888, 16.420, 29.467, 46.305, 67.170), nrow=5, ncol=3)
      model <- "without linear trend in shift correction"
    }
    N <- nrow(Z0)
    M00 <- crossprod(Z0)/N
    M11 <- crossprod(Z1)/N
    MKK <- crossprod(ZK)/N
    M01 <- crossprod(Z0, Z1)/N
    M0K <- crossprod(Z0, ZK)/N
    MK0 <- crossprod(ZK, Z0)/N
    M10 <- crossprod(Z1, Z0)/N
    M1K <- crossprod(Z1, ZK)/N
    MK1 <- crossprod(ZK, Z1)/N
    M11inv <- solve(M11)
    R0 <- Z0 - t(M01 %*% M11inv %*% t(Z1))
    RK <- ZK - t(MK1 %*% M11inv %*% t(Z1))
    S00 <- M00 - M01 %*% M11inv %*% M10
    S0K <- M0K - M01 %*% M11inv %*% M1K
    SK0 <- MK0 - MK1 %*% M11inv %*% M10
    SKK <- MKK - MK1 %*% M11inv %*% M1K
    Ctemp <- chol(SKK, pivot = TRUE)
    pivot <- attr(Ctemp, "pivot")
    oo <- order(pivot)
    C <- t(Ctemp[, oo])
    Cinv <- solve(C)
    S00inv <- solve(S00)
    valeigen <- eigen(Cinv %*% SK0 %*% S00inv %*% S0K %*% t(Cinv))
    lambda <- valeigen$values
    e <- valeigen$vector
    V <- t(Cinv) %*% e
    rownames(V) <- colnames(x)
    Vorg <- V
    V <- sapply(1:P, function(x) V[, x]/V[1, x])
    W <- S0K %*% V %*% solve(t(V) %*% SKK %*% V)
    PI <- S0K %*% solve(SKK)
    DELTA <- S00 - S0K %*% V %*% solve(t(V) %*% SKK %*% V) %*% 
        t(V) %*% SK0
    GAMMA <- M01 %*% M11inv - PI %*% MK1 %*% M11inv
    type <- "trace statistic"
    teststat <- as.matrix(rev(sapply(idx, function(x) N * sum(log(1 + lambda[(x + 1):P])))))
    colnames(teststat) <- "trace"
    if (arrsel > 5) {
      cval <- NULL
    } else {
      cval <- round(cvals[1:arrsel, ], 2)
      colnames(cval) <- c("10pct", "5pct", "1pct")
      rownames(cval) <- c(paste("r <= ", (arrsel - 1):1, " |", sep = ""), "r = 0  |")
    }
    temp1 <- NULL
    for (i in 1:(K - 1)) {
      temp <- paste(colnames(x), ".dl", i, sep = "")
      temp1 <- c(temp1, temp)
    }
    colnames(Z1) <- temp1
    colnames(ZK) <- paste(colnames(x), "l1", sep=".")
    colnames(Z0) <- paste(colnames(x), "d", sep=".")
    colnames(V) <- colnames(ZK)
    rownames(V) <- colnames(ZK)
    colnames(W) <- colnames(V)
    rownames(W) <- colnames(x)
    colnames(Vorg) <- colnames(V)
    rownames(Vorg) <- rownames(V)
    rownames(PI) <- colnames(x)
    colnames(PI) <- colnames(W)
    colnames(R0) <- paste("R0", colnames(Z0), sep = ".")
    colnames(RK) <- paste("RK", colnames(ZK), sep = ".")
    
    new("ca.jo", x = x, Z0 = Z0, Z1 = Z1, ZK = ZK, type = type, model = model, ecdet = "none", lag = K, P = arrsel, 
        season = season, dumvar = NULL, cval = cval, teststat = as.vector(teststat), 
        lambda = lambda, Vorg = Vorg, V = V, W = W, PI = PI, 
        DELTA = DELTA, GAMMA = GAMMA, R0 = R0, RK = RK, bp = tau.bp, 
        test.name = "Johansen-Procedure")
}
