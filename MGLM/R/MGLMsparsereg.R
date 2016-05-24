# Variable selection 
# Author: Zhang
#============================================================## 

#============================================================## 
## Set class
#============================================================##
setClass("MGLMsparsereg", representation(call = "function", data = "list", coefficients = "matrix", 
    logL = "numeric", BIC = "numeric", AIC = "numeric", Dof = "numeric", iter = "numeric", 
    maxlambda = "numeric", lambda = "numeric", distribution = "character", penalty = "character", 
    Beta = "numeric"))

##============================================================## 
## Sparse reg function 
##============================================================##

MGLMsparsereg <- function(formula, data, dist, lambda, penalty, weight, init, penidx, 
    maxiters = 150, ridgedelta, epsilon = 1e-05, regBeta = FALSE, overdisp) {
    call <- match.call()
    if (missing(data)) 
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "weight"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame(n = 1))
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "any")
    X <- model.matrix(mt, mf, contrasts)
    if (!missing(weight)) 
        weight <- weight[rowSums(Y) != 0]
    X <- as.matrix(X[rowSums(Y) != 0, ])
    Y <- as.matrix(Y[rowSums(Y) != 0, colSums(Y) != 0])
    d <- ncol(Y)
    m <- rowSums(Y)
    p <- ncol(X)
    N <- nrow(Y)
    if (missing(weight)) 
        weight <- rep(1, N)
    if (dist == "GDM" && d == 2) 
        stop("When d=2, GDM model is equivilant to DM model, please use dist='DM'.")
    if (penalty == "group") 
        penalty <- "group_row"
    if (!penalty %in% c("sweep", "group_row", "nuclear")) 
        stop("penalty type can only be sweep, group, or nuclear.")
    if (missing(penidx)) 
        penidx <- rep(TRUE, p)
    if (missing(init)) {
        if (dist == "MN") {
            init <- matrix(0, p, (d - 1))
        } else if (dist == "DM") {
            init <- matrix(0, p, d)
        } else if (dist == "GDM") {
            init <- matrix(0, p, 2 * (d - 1))
        } else if (dist == "NegMN") {
            if (regBeta) 
                init <- matrix(0, p, (d + 1)) else init <- matrix(0, p, d)
        }
    }
    if (dist == "NegMN" && regBeta == FALSE && missing(overdisp)) {
        est <- DMD.NegMN.fit(Y)
        overdisp <- est$estimate[d + 1]
    } else {
        overdisp <- NULL
    }
    if (missing(ridgedelta)) 
        ridgedelta <- 1/(p * d)
    est <- eval(call("MGLMsparsereg.fit", Y = Y, X = X, dist = dist, lambda = lambda, 
        penalty = penalty, weight = weight, init = init, penidx = penidx, maxiters = maxiters, 
        ridgedelta = ridgedelta, epsilon = epsilon, regBeta = regBeta, overdisp = overdisp))
    est$call <- match.call()
    est$data <- list(Y = Y, X = X)
    est$distribution <- ifelse(dist == "MN", "Multinomial", ifelse(dist == "DM", 
        "Dirichlet Multinomial", ifelse(dist == "GDM", "Generalized Dirichlet Multinomial", 
            "Negative Multinomial")))
    est$penalty <- ifelse(penalty == "group_row", "group", ifelse(penalty == "nuclear", 
        "nuclear", "sweep"))
    est$lambda <- lambda
    class(est) <- "MGLMsparsereg"
    return(est)
}

## ============================================================## 
## The MGLM sparse reg function 
## ============================================================##
MGLMsparsereg.fit <- function(Y, X, dist, lambda, penalty, weight, init, penidx, 
    maxiters = 150, ridgedelta, epsilon = 1e-05, regBeta = FALSE, overdisp) {
    
    d <- ncol(Y)
    m <- rowSums(Y)
    p <- ncol(X)
    N <- nrow(Y)
    
    if (dist == "GDM") {
        colOrder <- order(colSums(Y), decreasing = TRUE)
        Y <- Y[, colOrder]
        outOrder <- order(colOrder[-d])
        # Ys <- t( apply(apply( apply(Y,1,rev),2,cumsum),2,rev) )
    }
    
    beta_old <- init
    B <- init
    alpha <- 1
    alpha_iter <- list()
    objval <- Inf
    niter <- 1
    isdescent <- TRUE
    oriRidge <- ridgedelta
    while ((niter < 3) || (!stop)) {
        niter <- niter + 1
        beta_old <- B
        obj1 <- objval
        if (niter <= 2) 
            S <- init
        loss <- MGLM.loss(Y, X, S, dist, weight, regBeta, overdisp)
        loss.S <- loss[[1]]
        loss.D1S <- loss[[2]]
        for (l in 1:50) {
            A <- S - ridgedelta * loss.D1S
            B <- matrix(NA, nrow(B), ncol(B))
            B[!penidx, ] <- A[!penidx, ]
            if (d > 2) {
                Apen <- A[penidx, ]
            } else if (d == 2) 
                Apen <- matrix(A[penidx, ], , 1)
            pen <- matrix_threshold(X = Apen, lambda = ridgedelta * lambda, penalty = penalty)
            B[penidx, ] <- pen[[1]]
            penval <- pen[[2]]
            if (all(abs(B[penidx, ]) == 0)) {
                if (penalty == "sweep") {
                  maxlambda <- max(abs(Apen))/ridgedelta
                } else if (penalty == "group_row") {
                  maxlambda <- max(sqrt(rowSums(Apen^2)))/ridgedelta
                } else if (penalty == "nuclear") {
                  sing.vec <- svd(Apen)$d
                  maxlambda <- sing.vec[1]/ridgedelta
                }
                if (lambda > maxlambda) {
                  B[penidx, ] <- 0
                  penval <- 0
                }
            } else {
                maxlambda <- NULL
            }
            objloss <- MGLM.loss(Y, X, B, dist, weight, regBeta, overdisp)
            loss.S2 <- objloss[[1]]
            loss.D1S2 <- objloss[[2]]
            objval <- loss.S2 + penval
            BminusS <- B - S
            surval <- loss.S + sum(loss.D1S * BminusS) + norm(BminusS, type = "F")^2/2/ridgedelta + 
                penval
            if (!is.na(objval) && !is.na(surval)) 
                if (objval <= surval) 
                  break else ridgedelta <- ridgedelta/2
        }
        alpha_old <- alpha
        alpha <- (1 + sqrt(4 + alpha_old^2))/2
        if (!is.na(objval) & objval <= obj1) {
            stop <- abs(obj1 - objval) < epsilon * (abs(obj1) + 1)
            S <- B + (alpha_old - 1)/alpha * (B - beta_old)
        } else {
            objval <- obj1
            if (isdescent) {
                isdescent <- FALSE
                stop <- FALSE
                S <- B + (alpha_old - 1)/alpha * (beta_old - B)
                B <- beta_old
            } else {
                stop <- TRUE
            }
        }
        if (niter >= maxiters) 
            stop <- TRUE
    }
    if (d > 2) 
        Bpen <- B[penidx, ] else if (d == 2) 
        Bpen <- matrix(B[penidx, ], , 1)
    if (penalty == "sweep") {
        Dof <- sum(penidx == FALSE) * ncol(B) + sum(Bpen != 0)
    } else if (penalty == "group_row") {
        Dof <- sum(!penidx) * ncol(B) + sum(rowSums(Bpen^2) != 0) + sum((ncol(B) - 
            1) * rowSums(Bpen^2)/rowSums(Apen^2))
    } else if (penalty == "nuclear") {
        Aspectrum <- svd(A[penidx, ])$d
        if (sum(Aspectrum > ridgedelta * lambda) == 0) {
            Dof <- 0
        } else {
            Dof <- 0
            for (i in 1:sum(Aspectrum > ridgedelta * lambda)) {
                Dof <- Dof + 1 + 2 * sum(Aspectrum[i] * (Aspectrum[i] - ridgedelta * 
                  lambda)/(Aspectrum[i]^2 - Aspectrum[-i]^2)) + abs(p - d) * (1 - 
                  ridgedelta * lambda/Aspectrum[i])
            }
            Dof <- Dof + sum(!penidx) * d
        }
    }
    if (dist == "GDM") 
        Dof <- Dof/2 else Dof <- Dof
    logL <- -MGLM.loss(Y, X, B, dist, weight, regBeta, overdisp)[[1]]
    AIC <- -2 * logL + 2 * Dof
    BIC <- -2 * logL + log(N) * Dof
    
    if (dist == "GDM") {
        B <- B[, c(outOrder, outOrder + d - 1)]
    }
    
    return(list(coefficients = B, logL = logL, AIC = AIC, BIC = BIC, Dof = Dof, iter = niter, 
        maxlambda = maxlambda))
}

## ============================================================================##
## MGLM logss
## ============================================================================##

MGLM.loss <- function(Y, X, beta, dist, weight, regBeta = FALSE, Beta) {
    N <- nrow(Y)
    d <- ncol(Y)
    p <- ncol(X)
    m <- rowSums(Y)
    if (missing(weight)) 
        weight <- rep(1, N)
    if (dist == "MN") {
        P <- matrix(NA, N, d)
        P[, d] <- rep(1, N)
        P[, 1:(d - 1)] <- exp(X %*% beta)
        P <- P/rowSums(P)
        loss <- -sum(weight * dmn(Y, P))
        kr1 <- Y[, -d] - P[, 1:(d - 1)] * m
        if (d > 2) {
            lossD1 <- -colSums(kr(kr1, X, weight))
        } else if (d == 2) {
            lossD1 <- -colSums(kr1 * weight * X)
        }
        lossD1 <- matrix(lossD1, p, (d - 1))
    } else if (dist == "DM") {
        alpha <- exp(X %*% beta)
        loss <- -sum(weight * ddirm(Y, alpha))
        tmpvector <- digamma(rowSums(alpha) + m) - digamma(rowSums(alpha))
        tmpmatrix <- digamma(alpha + Y) - digamma(alpha)
        dalpha <- tmpmatrix - tmpvector
        lossD1 <- -rowSums(sapply(1:nrow(Y), function(i, A, B, w) return(w[i] * A[i, 
            ] %x% B[i, ]), alpha * dalpha, X, weight))
        lossD1 <- matrix(lossD1, p, d)
    } else if (dist == "GDM") {
        Ys <- t(apply(apply(apply(Y, 1, rev), 2, cumsum), 2, rev))
        alpha <- exp(X %*% beta)
        A <- alpha[, 1:(d - 1)]
        B <- alpha[, d:(2 * (d - 1))]
        loss <- -sum(weight * dgdirm(Y, A, B))
        dalpha1 <- digamma(A + Y[, -d]) - digamma(A) - digamma(A + B + Ys[, -d]) + 
            digamma(A + B)
        dalpha2 <- digamma(B + Ys[, -1]) - digamma(B) - digamma(A + B + Ys[, -d]) + 
            digamma(A + B)
        dalpha <- cbind(dalpha1, dalpha2)
        lossD1 <- -rowSums(sapply(1:nrow(Y), function(i, A, B, w) return(w[i] * A[i, 
            ] %x% B[i, ]), alpha * dalpha, X, weight))
        lossD1 <- matrix(lossD1, p, 2 * (d - 1))
    } else if (dist == "NegMN") {
        if (regBeta) {
            P <- matrix(NA, N, (d + 1))
            alpha <- exp(X %*% beta)
            Beta <- alpha[, (d + 1)]
            alpha_rowsums <- rowSums(alpha[, 1:d]) + 1
            P[, (d + 1)] <- 1/alpha_rowsums
            P[, 1:d] <- alpha[, 1:d] * P[, (d + 1)]
            loss <- -sum(weight * dneg(Y, alpha[, 1:d], Beta))
            deta <- matrix(0, nrow(Y), d + 1)
            deta[, 1:d] <- Y - alpha[, 1:d] * Beta - P[, 1:d] * (m - (alpha_rowsums - 
                1) * Beta)
            deta[, (d + 1)] <- Beta * (digamma(Beta + m) - digamma(Beta) + log(P[, 
                (d + 1)]))
            lossD1 <- -rowSums(sapply(1:nrow(Y), function(i, A, B, w) return(w[i] * 
                A[i, ] %x% B[i, ]), deta, X, weight))
            lossD1 <- matrix(lossD1, p, (d + 1))
        } else {
            P <- matrix(NA, N, (d + 1))
            alpha <- exp(X %*% beta)
            alpha_rowsums <- rowSums(alpha) + 1
            P[, (d + 1)] <- 1/alpha_rowsums
            P[, 1:d] <- alpha[, 1:d] * P[, (d + 1)]
            loss <- -sum(weight * dneg(Y, alpha, rep(Beta, N)))
            deta <- Y - (Beta + m) * P[, 1:d]
            lossD1 <- -rowSums(sapply(1:N, function(i, A, B, w) return(w[i] * A[i, 
                ] %x% B[i, ]), deta, X, weight))
            lossD1 <- matrix(lossD1, p, d)
        }
    }
    return(list(loss, lossD1))
}

## ============================================================================##
## matrix thresholding
## ============================================================================##

matrix_threshold <- function(X, lambda, penalty) {
    
    N <- nrow(X)
    d <- ncol(X)
    B <- matrix(0, N, d)
    if (penalty == "sweep") {
        B <- lsq_thresholding(X, lambda)
        penalty_value <- lambda * sum(abs(B))
    } else if (penalty == "group_row" || penalty == "group") {
        row_12norm <- sqrt(rowSums(X^2))
        vec <- 1 - lambda/row_12norm
        vec[vec < 0] <- 0
        B <- X * vec
        penalty_value <- lambda * sum(sqrt(rowSums(B^2)))
    } else if (penalty == "group_col") {
        row_12norm <- sqrt(colSums(X^2))
        B <- X * max(c(1 - lambda/row_12norm, 0))
        penalty_value <- lambda * sum(sqrt(colSums(B^2)))
    } else if (penalty == "nuclear") {
        decomp <- svt(X, lambda)
        U <- decomp[[1]]
        s <- as.vector(decomp[[2]])
        V <- decomp[[3]]
        if (length(s) == 0) 
            B <- matrix(0, N, d) else B <- U %*% (s * t(V))
        bs <- svd(B)$d
        penalty_value <- lambda * sum(bs)
    }
    
    return(list(B, penalty_value))
}

## ============================================================================##
## lsq_threshold We can only work on lasso for now
## ============================================================================##
lsq_thresholding <- function(b, lambda) {
    if (lambda < 0) 
        stop("penalty constant lambda should be nonnegative")
    
    B <- b
    B[abs(b) <= lambda] <- 0
    B[b > lambda] <- B[b > lambda] - lambda
    B[b < -lambda] <- B[b < -lambda] + lambda
    
    return(B)
}

## ============================================================================##
## svt
## ============================================================================##
svt <- function(b, lambda) {
    decomp <- svd(b)
    s <- decomp$d
    if (lambda > 0) {
        s <- lsq_thresholding(as.matrix(s), lambda)
        idx <- s > 0
        s <- s[idx]
        U <- decomp$u[, idx]
        V <- decomp$v[, idx]
    }
    return(list(U = U, s = s, V = V))
} 
