reReg <- function(formula, data, subset, method = "jsc", B = 200, contrasts = NULL) {
    Call <- match.call()
    ## obj <- eval(formula[[2]], data)
    mnames <- c("", "formula", "data", "subset")
    cnames <- names(Call)
    cnames <- cnames[match(mnames, cnames, 0)]
    mcall <- Call[cnames]
    mcall[[1]] <- as.name("model.frame")
    obj <- eval(mcall, parent.frame())
    df <- model.extract(obj, "response")
    if (!is.reSurv(df))
        stop("Response must be a reSurv object")
    formula[[2]] <- NULL
    if (formula == ~ 1) {
        df <- cbind(df, zero=0)
        method = "jsc"
    } else {
        df <- cbind(df, model.matrix(attr(obj, "terms"), obj, contrasts)[,-1])
        colnames(df)[-c(1:3)] <- colnames(model.matrix(attr(obj, "terms"), obj, contrasts))[-1]
    }
    df <- as.data.frame(df)
    df <- df[order(df$id, df$Time), ]
    id <- df$id
    delta <- df$event
    cluster <- unlist(lapply(split(id, id), function(x) 1:length(x)))
    clsz <- unlist(lapply(split(id, id), length))
    mt <- unlist(lapply(split(cluster, id), length)) - 1
    Y <- rep(df$Time[cumsum(clsz)], clsz)
    T <- df$Time
    X <- as.matrix(df[,-c(1:3)])
    row.names(X) <- NULL
    out <- NULL
    out <- reReg.fit(X = X, Y = Y, T = T, id = id, cluster = cluster, delta = delta, B = B, method = method)
    ## Bootstrap method for method != jsc
    if (B > 0 & method != "jsc") {
        p <- ncol(X)
        n <- length(unique(id))
        bout <- matrix(0, ncol = 2 * p, nrow = B)
        for (i in 1:B) {
            sampled.id <- sample(unique(id), n, TRUE)
            ind <- unlist(sapply(sampled.id, function(x) which(id == x)))
            new.id <- rep(1:n, clsz[sampled.id])
            tmp <- reReg.fit(X = as.matrix(X[ind,]), Y = Y[ind], T = T[ind], id = new.id,
                             cluster = unlist(lapply(split(new.id, new.id), function(x) 1:length(x))),
                             delta = delta[ind], B = 0, method = method)
            bout[i,] <- c(tmp$alpha, tmp$beta)
        }
        out$va <- var(bout[,1:p])
        out$vb <- var(bout[,1:p + p])
        out$aSE <- sqrt(diag(out$va))
        out$bSE <- sqrt(diag(out$vb))
    }
    ## haz/rate estimations
    alpha <- out$alpha
    beta <- out$beta
    Ya <- log(Y) + X %*% alpha
    Ta <- log(T) + X %*% alpha
    ord <- order(T)
    Ti <- T[ord]
    dummy <- 1:length(T)
    repeats <- table(Ti)
    lambda <- npMLE(unique(log(Ti)), Ta, Ya)
    lambda <- rep(lambda, repeats)
    lambda[dummy[ord]] <- lambda / max(lambda)
    ltmp <- npMLE(Ya[cluster == 1], Ta, Ya)
    zHat <- as.numeric(mt * max(ltmp) / ltmp)
    zHat <- ifelse(zHat %in% c("Inf", "NA", "NaN"), 0, zHat)
    Yb <- log(Y) + X %*% beta
    Yb <- Yb[cluster == 1]
    haz <- sapply(unique(Ti), function(x) baseHaz(x, exp(Yb), zHat, delta[cluster == 1]))
    haz <- rep(haz, repeats)
    haz[dummy[ord]] <- haz
    out$lambda <- lambda
    out$haz <- haz
    out$zHat <- zHat
    out$df <- df
    out$Call <- Call
    out$method <- method
    names(out$alpha) <- names(out$beta) <- names(df)[-c(1:3)]
    class(out) <- "reReg"
    out
}

npMLE <- function(t, tij, yi, weights = NULL) {
    if (is.null(weights))
        weights <- rep(1, length(yi))
    ttmp <- tij[tij != yi]
    ord <- order(ttmp)
    sl <- unique(ttmp[ord])
    l <- ifelse(min(t) < max(sl), which(sl > min(t))[1], length(sl))
    ## res <- vector("double", 1)
    ## tmp <- sl[l:length(sl)]
    tmp <- sl[l:length(sl)]
    tmp <- rev(tmp)
    tij <- rev(tij)
    yi <- rev(yi)
    ## yi <- ifelse(is.infinite(yi), max(yi[!is.infinite(yi)]), yi)
    ## tij <- ifelse(is.infinite(tij), max(tij[!is.infinite(tij)]), tij)
    res <- vector("double", length(tmp)) + 1
    res <- .C("plLambda", as.double(tmp), as.double(tij), as.double(yi), as.double(weights), 
              as.integer(length(tmp)), as.integer(length(yi)),
              out = as.double(res), PACKAGE = "reReg")$out
    out <- rev(res)[sapply(t, function(x) which(rev(tmp) >= x)[1])]
    out <- ifelse(is.na(out), 0, out)
    out <- exp(-out)
}


baseHaz <- function(t, Y, zHat, delta, weights  = NULL) {
    if (is.null(weights)) 
        weights <- rep(1, length(Y))
    ind <- which(delta == 1 & Y <= t)
    temp2 <- tmp <- weights[order(Y)]
    ## temp2 <- c(tmp[1], diff(cumsum(tmp)))
    ## temp2[order(Y)] <- temp2
    temp2[order(Y)] <- tmp
    if (length(ind) > 0) {
        out <- sapply(ind, function(x) temp2[x] / sum(zHat * weights * (Y >= Y[x])))
    }
    if (length(ind) == 0)
        out <- 0
    sum(out)
}

alphaEq <- function(alpha, X, Y, T, cluster, mt, esteq = "jsc", weights = NULL) {
    n <- sum(cluster == 1)
    if (is.null(weights))
        weights <- rep(1, n)
    p <- ncol(X)
    ## Ystar <- Y * exp(X %*% alpha)
    ## Tstar <- T * exp(X %*% alpha)
    Ystar <- log(Y) + X %*% alpha
    Tstar <- log(T) + X %*% alpha
    Lambda <- npMLE(Ystar[which(cluster == 1)], Tstar, Ystar,
                    weights = rep(weights, diff(c(which(cluster ==1), length(cluster)+1))))
    ## Lambda <- npMLE(Ystar[which(cluster == 1)], log(T), log(Y),
    ## weights = rep(weights, diff(c(which(cluster ==1), length(cluster)+1))))
    res <- vector("double", p * length(weights) %/% n)
    if (esteq == "jsc") 
        res <- .C("alphaEq1", as.double(X[which(cluster == 1), ]), as.double(Lambda),
                  as.double(weights), as.integer(mt), as.integer(n), as.integer(p),
                  as.integer(length(weights) %/% n), 
                  out = as.double(res), PACKAGE = "reReg")$out
    if (esteq == "M2")
        res <- .C("alphaEq2", as.double(X[which(cluster == 1),]), as.double(Lambda),
                  as.integer(mt), as.integer(n), as.integer(p), out = as.double(res),
                  PACKAGE = "reReg")$out
    res / rep(n * unlist(lapply(split(weights, rep(1:(length(weights) %/% n), each = n)), sum)), each = p)
    ## res / n ^ 2
}

betaEq <- function(X, Y, T, cluster, mt, delta, zHat = NULL, alpha, beta, weights = NULL) {
    p <- ncol(X)
    n <- sum(cluster == 1)
    if (is.null(weights))
        weights <- rep(1, n)
    if (is.null(zHat)) {
        Ystar <- log(Y) + X %*% alpha
        Tstar <- log(T) + X %*% alpha
        lambda <- npMLE(Ystar[which(cluster == 1)], Tstar, Ystar,
                        weights = rep(weights, diff(c(which(cluster ==1), length(cluster)+1))))
        zHat <- as.numeric(weights * mt / lambda)
        zHat <- ifelse(zHat %in% c("Inf", "NA", "NaN"), 0, zHat)
    }
    Y <- log(Y) + X %*% beta
    Y <- Y[which(cluster == 1)]
    X <- X[which(cluster == 1), ]
    delta <- delta[which(cluster == 1)]
    res <- vector("double", p * length(weights) %/% n)
    res <- .C("betaEst", as.double(Y), as.double(X), as.double(delta), as.double(zHat),
              as.double(weights), as.integer(n), as.integer(p),
              as.integer(length(weights) %/% n), 
              out = as.double(res), PACKAGE = "reReg")$out
    res / n
}

ghoshU2 <- function(alpha, beta, T, Y, X, cl) {
    ## dim(X) = n by p, dim(Y) = n by 1, dim(T) > n by 1
    d <- max(X %*% (alpha - beta), 0)
    TT <- log(T) - rep(X %*% alpha, cl)
    TY <- log(Y) - X %*% beta - d
    p <- ncol(X)
    .C("ghosh", as.double(TT), as.double(TY), as.double(X), as.integer(cl),
       as.integer(c(0, cumsum(cl)[-length(cl)])),
       as.integer(nrow(X)), as.integer(p), 
       out = as.double(double(p)), PACKAGE = "reReg")$out
}

reReg.fit <- function(X, Y, T, id, cluster, delta, method = "jsc", B = 200) {
    n <- sum(cluster == 1)
    mt <- unlist(lapply(split(cluster, id), length)) - 1
    p <- ncol(X)
    ## Point estiamtion
    alpha <- beta <- gamma <- rep(0, p)
    if (method == "GL") {
        outB <- aftsrr(Surv(Y, delta) ~ X, subset = cluster == 1, B = 0,
                       rankWeights = "logrank", method = "nonsm")
        ## outB <- faft::faft(Y[cluster == 1], delta[cluster == 1], X[cluster == 1,])
        ## outA <- dfsane(alpha, ghoshU2, beta = outB$beta, T = T[T != Y],
        ##                Y = Y[cluster == 1 & T != Y],
        ##                X = as.matrix(X[cluster == 1 & T != Y, ]),
        ##                cl = unlist(lapply(split(id[T != Y], id[T != Y]), length)),
        ##                alertConvergence = FALSE, quiet = TRUE, 
        ##                control = list(NM = FALSE, M = 100, noimp = 50, trace = FALSE))
        outA <- dfsane(alpha, ghoshU2, beta = outB$beta, T = ifelse(T == Y, 1e5, T),
                       Y = Y[cluster == 1],
                       X = as.matrix(X[cluster == 1, ]),
                       cl = unlist(lapply(split(id, id), length)),
                       alertConvergence = FALSE, quiet = TRUE, 
                       control = list(NM = FALSE, M = 100, noimp = 50, trace = FALSE))
        outB$par <- -1 * outB$beta
        outA$par <- -1 * outA$par
        zHat <- NA
    }
    if (method == "jsc") {
        outA <- dfsane(alpha, alphaEq, X = X, Y = Y, T = T, cluster = cluster, mt = mt,
                       esteq = method, weights = NULL, alertConvergence = FALSE, quiet = TRUE, 
                       control = list(NM = FALSE, M = 100, noimp = 50, trace = FALSE))
        alpha <- outA$par
        Ystar <- log(Y) + X %*% alpha
        Tstar <- log(T) + X %*% alpha
        lambda <- npMLE(Ystar[cluster == 1], Tstar, Ystar)
        ## zHat <- as.numeric(mt / lambda)
        zHat <- as.numeric(mt * npMLE(log(max(Y)), Tstar, Ystar) / lambda)
        zHat <- ifelse(zHat %in% c("Inf", "NA", "NaN"), 0, zHat)
        outB <- dfsane(beta, betaEq, X = X, Y = Y, T = T, cluster = cluster, delta = delta, mt = mt,
                       alpha = outA$par, zHat = zHat, weights = NULL, alertConvergence = FALSE, quiet = TRUE, 
                       control = list(NM = FALSE, M = 100, noimp = 50, trace = FALSE))
    }
    if (method == "HW") {
        gamma <- c(0, gamma)
        ## outA <- BBsolve(gamma,
        X <- cbind(1, X[cluster == 1,])
        outA <- dfsane(gamma, HWeq, X = X, Y = Y, T = T, cluster = cluster, mt = mt, alertConvergence = FALSE,
                       quiet = TRUE, control = list(NM = FALSE, M = 100, noimp = 50, trace = FALSE))
        alpha <- outA$par <- outA$par[-1]
        lambda <- npMLE(Y[cluster == 1], T, Y)
        zHat <- as.numeric(mt * npMLE(max(Y), T, Y) / (lambda * exp(X[, -1] %*% alpha)))
        zHat <- ifelse(zHat %in% c("Inf", "NA", "NaN"), 0, zHat)
        outB <- dfsane(beta, HWeq2, X = X[, -1], Y = Y[cluster == 1], delta = delta[cluster == 1], zHat = zHat,
                       alertConvergence = FALSE, quiet = TRUE,
                       control = list(NM = FALSE, M = 100, noimp = 50, trace = FALSE))
    }
    ## Variance estimation
    aSE <- bSE <- da <- va <- db <- vb <- NA
    if (B > 0 & method == "jsc") {
        ## var
        E <- matrix(rexp(n * B), nrow = n)
        Z <- matrix(rnorm(p * B), nrow = p)
        ua <- matrix(apply(Z, 2, function(x) alphaEq(outA$par + n ^ (-0.5) * x, X, Y, T, cluster, mt)),
                     nrow = p)
        da <- t(apply(ua, 1, function(x) lm(n ^ (0.5) * x ~ t(Z))$coef[-1]))
        ## ua2 <- matrix(alphaEq(outA$par, X, Y, T, id, mt, weights = as.double(E)), nrow = p)
        ua2 <- apply(E, 2, function(x) alphaEq(outA$par, X, Y, T, cluster, mt, weights = x))
        va <- var(t(matrix(ua2, nrow = p)))
        if (qr(da)$rank == p)
            aVar <- solve(da) %*% va %*% t(solve(da))
        if (qr(da)$rank != p)
            aVar <- ginv(da) %*% va %*% t(ginv(da))
        aSE <- sqrt(diag(aVar))
        ub <- matrix(apply(Z, 2, function(x)
            betaEq(X, Y, T, cluster, mt, delta, zHat, outA$par, outB$par + n ^ (-0.5) * x)), nrow = p)
        db <- t(apply(ub, 1, function(x) lm(n ^ (0.5) * x ~ t(Z))$coef[-1]))
        ## ub2 <- matrix(betaEq(X, Y, T, id, mt, delta, NULL, outB$par, weights = as.double(E)), nrow = p)
        ub2 <- apply(E, 2, function(x) betaEq(X, Y, T, cluster, mt, delta, NULL, outA$par, outB$par, weights = x))
        vb <- var(t(matrix(ub2, nrow = p)))
        if (qr(db)$rank == p)
            bVar <- solve(db) %*% vb %*% t(solve(db))
        if (qr(db)$rank != p)
            bVar <- ginv(db) %*% vb %*% t(ginv(db))
        bSE <- sqrt(diag(bVar))    
    }
    out <- list(alpha = outA$par, aSE = aSE, aconv = outA$convergence,
                beta = outB$par, bSE = bSE, bconv = outB$convergence,
                da = da, va = va, db = db, vb = vb, muZ = mean(zHat))
    out
    ## c(outA$par, aSE, outA$convergence, outB$par, bSE, outB$convergence)
}


##########################################################################################
## Paper 2: More general models
##########################################################################################

sarm <- function(X, Y, T, id, cluster, method, B = 200) {
    n <- sum(cluster == 1)
    mt <- unlist(lapply(split(cluster, id), length)) - 1
    p <- ncol(X)
    alpha <- beta <- gamma <- rep(0, p)
    muZ <- NULL
    if (method == "M1") {
        gamma <- c(0, gamma)
        out <- BBsolve(gamma, coefEq, alpha = alpha, X = X, Y = Y, T = T,
                       cluster = cluster, mt = mt, weights = NULL,
                       quiet = TRUE, control = list(M = c(1, 10)))
        muZ <- out$par[1]
        alpha <- out$par[-1]
    }
    if (method == "M3") {
        alpha <- BBsolve(alpha, M1eq, X = X, Y = Y, T = T, cluster = cluster, weights = NULL,
                         quiet = TRUE, control = list(M = c(1, 10)))$par
        gamma <- c(0, gamma)
        if (alpha %*% alpha > 100) {
            beta <- c(0, alpha)
        }
        else {
            out <- BBsolve(gamma, coefEq, alpha = alpha, X = X, Y = Y, T = T,
                           cluster = cluster, mt = mt, weights = NULL,
                           quiet = TRUE, control = list(M = c(1, 10)))
            muZ <- out$par[1]
            beta <- out$par[-1] + alpha
        }
    }
    if (method == "M2") {
        ## gamma <- c(0, gamma)
        ## out <- BBsolve(gamma, coefEq, alpha = NULL, X = X, Y = Y, T = T,
        ##                cluster = cluster, mt = mt, weights = NULL,
        ##                quiet = TRUE, control = list(M = c(1, 10)))
        ## muZ <- out$par[1]
        ## alpha <- out$par[-1]
        ## out <- dfsane(alpha, M1eq, X = X, Y = Y, T = T,
        ##                 cluster = cluster, weights = NULL,
        ##                 control = list(NM = TRUE, M = 1, noimp = 50, trace = FALSE))
        out <- BBsolve(alpha, M1eq, X = X, Y = Y, T = T, cluster = cluster, weights = NULL,
                       quiet = TRUE, control = list(M = c(1, 10)))
        alpha <- out$par
    }
    list(alpha = alpha, beta = beta, muZ = muZ)
}

HWeq <-function(gamma, X, Y, T, cluster, mt) {
    n <- sum(cluster == 1)
    Lambda <- npMLE(Y[cluster == 1], T, Y)
    res <- vector("double", length(gamma))
    p <- ncol(X)
    .C("sarm1", as.double(X), as.double(Lambda), as.double(rep(1, n)),
       as.double(gamma), as.integer(mt), as.integer(n), as.integer(p), as.integer(1),
       out = as.double(rep(0, p)), PACKAGE = "reReg")$out                            
}


HWeq2 <-function(beta, X, Y, delta, zHat) {
    n <- nrow(X)
    p <- ncol(X)
    res <- vector("double", p)
    res <- .C("HWb", as.double(Y), as.double(X), as.double(delta), as.double(zHat),
              as.double(X %*% beta), as.integer(n), as.integer(p), as.integer(1),
              out = as.double(res), PACKAGE = "reReg")$out
    res / n
}

coefEq <- function(alpha, gamma, X, Y, T, cluster, mt, weights = NULL) {
    n <- sum(cluster == 1)
    if (is.null(weights))
        weights <- rep(1, n)
    if (is.null(alpha))
        alpha <- -1 * gamma[-1]
    Ytmp <- log(Y) + X %*% alpha
    Ttmp <- log(T) + X %*% alpha    
    Lambda <- npMLE(Ytmp[cluster == 1], Ttmp, Ytmp,
                    weights = rep(weights, diff(c(which(cluster ==1), length(cluster)+1))))
    Lambda <- Lambda / max(Lambda)
    X <- X[cluster == 1,]
    p <- ncol(X)
    res <- vector("double", (p + 1) * length(weights) %/% n)
    res <- .C("sarm1", as.double(cbind(1, X)), as.double(Lambda), as.double(weights),
              as.double(gamma), as.integer(mt), as.integer(n), as.integer(p+1),
              as.integer(length(weights) %/% n),
              out = as.double(res), PACKAGE = "reReg")$out
    res / n    
}

## SARM 3
## M1eq <- function(alpha, X, Y, T, cluster, weights = NULL) {
##     n <- sum(cluster == 1)
##     ## mi <- unlist(lapply(split(cluster, id), length)) - 1
##     if (is.null(weights))
##         weights <- rep(1, n)
##     p <- ncol(X)
##     Ytmp <- log(Y) + X %*% alpha
##     Ttmp <- log(T) + X %*% alpha
##     ## Ytmp <- Y * exp(X %*% alpha)
##     ## Ttmp <- T * exp(X %*% alpha)
##     lambda <- npMLE(Ytmp[which(cluster == 1)], Ttmp, Ytmp)
##     ind <- which(T != Y)
##     Ttmp <- Ttmp[ind]
##     Ytmp <- Ytmp[ind]
##     X <- X[ind,]
##     mt <- c(diff(which(cluster == 1)), length(cluster) - rev(which(cluster==1))[1] + 1) - 1
##     lambda <- rep(lambda, mt)
##     mt <- rep(mt, mt)
##     res <- vector("double", p * length(weights) %/% n)
##     res <- .C("sarm2", as.double(X), as.double(Ttmp), as.double(Ytmp), as.double(weights),
##               as.double(lambda), as.double(mt), 
##               as.integer(length(Ttmp)), as.integer(p), as.integer(length(weights) %/% n),
##               out = as.double(res), PACKAGE = "reReg")$out
##     res
## }

## Martingal approach
M1eq <- function(alpha, X, Y, T, cluster, weights = NULL) {
    n <- sum(cluster == 1)
    ## mi <- unlist(lapply(split(cluster, id), length)) - 1
    if (is.null(weights))
        weights <- rep(1, n)
    p <- ncol(X)
    Ytmp <- log(Y) + X %*% alpha
    Ttmp <- log(T) + X %*% alpha
    ## Ytmp <- Y * exp(X %*% alpha)
    ## Ttmp <- T * exp(X %*% alpha)
    ind <- which(T != Y)
    Ttmp <- Ttmp[ind]
    Ytmp <- Ytmp[ind]
    X <- X[ind,]
    res <- vector("double", p * length(weights) %/% n)
    res <- .C("sarm2", as.double(X), as.double(Ttmp), as.double(Ytmp), as.double(weights), 
              as.integer(length(Ttmp)), as.integer(p), as.integer(length(weights) %/% n),
              out = as.double(res), PACKAGE = "reReg")$out
    res
}

## Method 1: Z\lambda(t)e^xa  CY's 2004 JASA
## Method 2: Z\lambda(te^xa)
## Method 3: Z\lambda(te^xa)e^xb
## Method 4: Z\lambda(te^xa)e^xa

##################################
## Newton's method
## still underconstruction
##################################

reReg2 <- function(X, Y, T, id, cluster, delta, method, B = 100, a1, b1, maxit = 100) {
    n <- sum(id == 1)
    mt <- unlist(lapply(split(id, cluster), length)) - 1
    p <- ncol(X)
    alpha <- beta <- rep(0, p)
    E <- matrix(rexp(n * B), nrow = n)
    Z <- matrix(rnorm(p * B), nrow = p)
    s <- 0
    for (i in 1:maxit) {
        uStar <- apply(Z, 2, function(x) alphaEq(a1 + n ^ (-0.5) * x, X, Y, T, id, mt))
        numDev <- t(apply(uStar, 1, function(x) lm(n ^ (0.5) * x ~ t(Z))$coef[-1]))
        ## score function
        sf <- solve(numDev) %*% alphaEq(a1, X, Y, T, id, mt)
        for (j in c(100, 0:10)) {
            a2 <- a1 - 2 ^ -j * sf
            u <- alphaEq(a2, X, Y, T, id, mt)
            ## uStar2 <- apply(E, 2, function(x) alphaEq(a2, X, Y, T, id, mt, weights = x))
            uStar2 <- matrix(alphaEq(a2, X, Y, T, id, mt, weights = as.double(E)), nrow = p)
            v <- var(t(uStar2))
            omega <- t(u) %*% solve(v) %*% u
            if (s >= omega) {
                s <- 0
                break
            }
            else
                s <- omega
        }
        if (!all(abs(a2 - a1) > 0.001))
            break
        else
            a1 <- a2
    }
    print(i)
    aVar <- sqrt(diag(solve(numDev) %*% v %*% t(solve(numDev))))
    c(a2, aVar)
}

