# This is package RAD 

"dlbetabinom" <-
function (nTot, n, alpha, alphaStar) 
{
    if (any(alphaStar == 0)) 
        return(NaN)
    else return(lchoose(nTot, n) + lgamma(n + alpha) + lgamma(nTot - 
        n + alphaStar) - lgamma(alpha + nTot + alphaStar) - (lgamma(alpha) + 
        lgamma(alphaStar) - lgamma(alpha + alphaStar)))
}


"dlDirichletMultinomial" <-
function (nID, pi, theta) 
{
    theta.j <- theta * pi
    nTot <- sum(nID)
    d <- lgamma(nTot + 1) - sum(lgamma(nID + 1)) + lgamma(theta) - 
        lgamma(nTot + theta) + sum(lgamma(nID + theta.j)) - sum(lgamma(theta.j))
    return(d)
}


"estimate.dirichletMultinomial" <-
function (n, ID, X, logij, est.var = FALSE, calc.resid = FALSE, 
    trace = TRUE) 
{
    print("Estimating parameters for the Dirichlet multinomial model")
    my.fun <- function(x) {
        -sumLogl(tau = x[-1], theta = x[1], nu = NULL, n, X, 
            logij, ID, "DMn")
    }
    inits <- rep(1e-04, ncol(X) + 1)
    cat("Ite:     -logl    : ", paste("theta   :  ", paste(colnames(X), 
        collapse = " :   ")), "\n")
    lower <- c(0, rep(-Inf, ncol(X)))
    upper <- rep(Inf, ncol(X) + 1)
    esti <- nlminb(start = inits, my.fun, control = list(trace = TRUE), 
        lower = lower, upper = upper)
    print(esti$message)
    names(esti$par) <- c("theta", colnames(X))
    AIC <- 2 * length(esti$par) + 2 * esti$objective
    var <- 0
    if (est.var) {
        print("Calculating the variance of the estiamtes")
        var <- solve(nH2(pt = esti$par, fun = my.fun))
        colnames(var) <- rownames(var) <- names(esti$par)
    }
    theta <- esti$par[1]
    Qres <- 0
    my.fit <- NULL
    if (calc.resid) {
        print("Calculating expectations and residuals")
        pij <- exp(-(X %*% esti$par[-1]) * logij)
        q.res <- mu <- rep(0, nrow(X))
        for (ii in unique(ID)) {
            pij[ID == ii] <- pij[ID == ii]/sum(pij[ID == ii])
            N <- sum(n[ID == ii])
            S <- length(n[ID == ii])
            mu[ID == ii] <- pij[ID == ii] * N
            alpha <- pij[ID == ii] * theta
            alpha.star <- sum(alpha) - alpha
            alpha.star <- ifelse(alpha.star < 0, 0, alpha.star)
            for (kk in 1:S) q.res[ID == ii][kk] <- sum(exp(dlbetabinom(N, 
                0:(n[ID == ii][kk] - 1), alpha[kk], alpha.star[kk])), 
                exp(dlbetabinom(N, n[ID == ii][kk], alpha[kk], 
                  alpha.star[kk]))/2)
        }
        Qres <- qnorm(q.res)
        my.fit <- cbind(mu, pij)
        colnames(my.fit) <- c("mu", "prob")
    }
    return(list(coef = esti$par, vcov = var, logl = -esti$objective, 
        AIC = AIC, fitted = my.fit, residuals = Qres))
}


"estimate.modifiedDMn" <-
function (n, ID, X, logij, est.var = FALSE, calc.resid = FALSE, 
    trace = TRUE) 
{
    print("Using the multinomial model to get initial values for estimation")
    esti <- estimate.multinomial(n, ID, X, logij, trace = trace)
    print("Estimating parameters for the modified Dirichlet multinomial model")
    my.fun <- function(x) {
        -sumLogl(tau = x[-(1:2)], theta = x[1], nu = x[2], n, 
            X, logij, ID, "MDMn")
    }
    inits <- c(1e-04, 1e-04, esti$coef)
    cat("Ite:     -logl    : ", paste("theta  :  ", "nu      :  ", 
        paste(colnames(X), collapse = "  :   ")), "\n")
    lower <- c(0.001, 0.001, rep(-Inf, ncol(X)))
    upper <- c(Inf, Inf, rep(Inf, ncol(X)))
    esti <- nlminb(start = inits, my.fun, control = list(trace = trace, 
        eval.max = 1000, iter.max = 400), lower = lower, upper = upper)
    print(esti$message)
    names(esti$par) <- c("theta", "nu", colnames(X))
    AIC <- 2 * length(esti$par) + 2 * esti$objective
    var <- 0
    if (est.var) {
        print("Calculating the variance of the estiamtes")
        var <- solve(nH2(pt = esti$par, fun = my.fun))
        colnames(var) <- rownames(var) <- names(esti$par)
    }
    theta <- esti$par[1]
    nu <- esti$par[2]
    Qres <- 0
    my.fit <- NULL
    if (calc.resid) {
        print("Calculating expectations and residuals")
        pij <- exp(-(X %*% esti$par[-(1:2)]) * logij)
        q.res <- mu <- rep(0, nrow(X))
        for (ii in unique(ID)) {
            pij[ID == ii] <- pij[ID == ii]/sum(pij[ID == ii])
            N <- sum(n[ID == ii])
            S <- length(n[ID == ii])
            mu[ID == ii] <- pij[ID == ii] * N
            alpha <- pij[ID == ii] * theta
            alpha.star <- sum(alpha) - alpha
            q.res[ID == ii] <- .Call("estBetaBinomialQuantile", 
                N, n[ID == ii], alpha, alpha.star, 1, -nu, PACKAGE = "RAD")
        }
        Qres <- q.res
        my.fit <- data.frame(nij = mu, pij = pij)
    }
    return(list(coef = esti$par, vcov = var, logl = -esti$objective, 
        AIC = AIC, fitted = my.fit, residuals = Qres, convergence = esti$convergence))
}


"estimate.multinomial" <-
function (n, ID, X, logij, est.var = FALSE, calc.resid = FALSE, 
    trace = TRUE) 
{
    print("Estimating parameters for the multinomial model")
    my.fun <- function(x) {
        -sumLogl(tau = x, theta = NULL, nu = NULL, n, X, logij, 
            ID, "Mn")
    }
    inits <- rep(1e-04, ncol(X))
    cat("Ite:     -logl    : ", paste(colnames(X), collapse = " :   "), 
        "\n")
    esti <- nlminb(start = inits, my.fun, control = list(trace = trace, 
        iter.max = 400))
    print(esti$message)
    names(esti$par) <- colnames(X)
    AIC <- 2 * length(esti$par) + 2 * esti$objective
    var <- 0
    if (est.var) {
        print("Calculating the variance of the estiamtes")
        var <- solve(nH2(pt = esti$par, fun = my.fun))
        colnames(var) <- rownames(var) <- names(esti$par)
    }
    Qres <- 0
    my.fit <- NULL
    if (calc.resid) {
        print("Calculating expectations and residuals")
        pij <- exp(-(X %*% esti$par) * logij)
        q.res <- mu <- rep(0, nrow(X))
        for (ii in unique(ID)) {
            pij[ID == ii] <- pij[ID == ii]/sum(pij[ID == ii])
            N <- sum(n[ID == ii])
            S <- length(n[ID == ii])
            mu[ID == ii] <- pij[ID == ii] * N
            for (kk in 1:S) q.res[ID == ii][kk] <- sum(dbinom(0:(n[ID == 
                ii][kk] - 1), N, pij[ID == ii][kk]), dbinom(n[ID == 
                ii][kk], N, pij[ID == ii][kk])/2)
        }
        Qres <- qnorm(q.res)
        my.fit <- cbind(mu, pij)
        colnames(my.fit) <- c("mu", "prob")
    }
    return(list(coef = esti$par, vcov = var, logl = -esti$objective, 
        AIC = AIC, fitted = my.fit, residuals = Qres))
}


"MDMnMod" <-
function (MDMn.form, data, ID, dist = "MDMn", scale.covar = FALSE, 
    est.var = TRUE, calc.resid = TRUE, trace = TRUE) 
{
    temp <- model.frame(MDMn.form, as.data.frame(data))
    n <- model.response(temp)
    names(n) <- NULL
    X <- model.matrix(MDMn.form, data)
    mean.X <- sd.X <- NA
    if (scale.covar & ncol(X) > 1) {
        X.t <- data.frame(X)
        mean.X <- apply(as.data.frame(X.t[, 2:ncol(X.t)]), 2, 
            mean)
        names(mean.X) <- names(X.t)[2:ncol(X.t)]
        sd.X <- apply(X.t[, 2:ncol(X.t)], 2, sd)
        X[, 2:ncol(X)] <- scale(X[, 2:ncol(X)])
    }
    logij <- rep(0, length(ID))
    for (ii in unique(ID)) logij[ID == ii] <- log(1:length(ID[ID == 
        ii]))
    if (dist == "Mn") {
        esti <- estimate.multinomial(n, ID, X, logij, est.var, 
            calc.resid, trace)
    }
    if (dist == "DMn") {
        esti <- estimate.dirichletMultinomial(n, ID, X, logij, 
            est.var, calc.resid, trace)
    }
    if (dist == "MDMn") {
        esti <- estimate.modifiedDMn(n, ID, X, logij, est.var, 
            calc.resid, trace)
    }
    esti$mean.X <- mean.X
    esti$sd.X <- sd.X
    esti$formula <- MDMn.form
    class(esti) <- "MDMnMod"
    esti
}


"NBLogl" <-
function (tau, od, X, y, offset) 
{
    lp <- X %*% tau + offset
    lam <- exp(lp)
    lnb <- dnbinom(y, size = od, mu = lam, log = TRUE)
    return(sum(lnb))
}


"nd2" <-
function (x0, f, m = NULL, D.accur = 4, ...) 
{
    D.n <- length(x0)
    if (is.null(m)) {
        D.f0 <- f(x0, ...)
        m <- length(D.f0)
    }
    if (D.accur == 2) {
        D.w <- tcrossprod(rep(1, m), c(-1/2, 1/2))
        D.co <- c(-1, 1)
    }
    else {
        D.w <- tcrossprod(rep(1, m), c(1/12, -2/3, 2/3, -1/12))
        D.co <- c(-2, -1, 1, 2)
    }
    D.n.c <- length(D.co)
    macheps <- .Machine$double.eps
    D.h <- macheps^(1/3) * abs(x0)
    D.deriv <- matrix(NA, nrow = m, ncol = D.n)
    for (ii in 1:D.n) {
        D.temp.f <- matrix(0, m, D.n.c)
        for (jj in 1:D.n.c) {
            D.xd <- x0 + D.h[ii] * D.co[jj] * (1:D.n == ii)
            D.temp.f[, jj] <- f(D.xd, ...)
        }
        D.deriv[, ii] <- rowSums(D.w * D.temp.f)/D.h[ii]
    }
    return(as.double(D.deriv))
}


"negBinMod" <-
function (NB.form, data, est.var = TRUE, scale.covar = FALSE, 
    trace = TRUE) 
{
    my.fun <- function(x) {
        -NBLogl(x[-1], x[1], X, y, offset)
    }
    environment(NB.form) <- environment()
    environment(data) <- environment(NB.form)
    temp <- model.frame(NB.form, as.data.frame(data))
    y <- model.response(temp)
    names(y) <- NULL
    X <- model.matrix(NB.form, data)
    mean.X <- sd.X <- NA
    if (scale.covar) {
        X.t <- data.frame(X)
        mean.X <- apply(X.t[, 2:ncol(X.t)], 2, mean)
        sd.X <- apply(X.t[, 2:ncol(X.t)], 2, sd)
        X[, 2:ncol(X)] <- scale(X[, 2:ncol(X)])
    }
    offset <- model.offset(temp)
    if (is.null(offset)) 
        offset <- rep(0, nrow(temp))
    print("Estimating parameters")
    print("Finding initial values using profile quasi - likelihood")
    fm.nb <- glm.nb(NB.form, data = data, init.theta = 1)
    inits <- c(fm.nb$theta, fm.nb$coef)
    print("Final estimation")
    cat("Ite:     -logl    :  Disp    ", paste(colnames(X), collapse = " :   "), 
        "\n")
    esti <- nlminb(start = inits, my.fun, control = list(trace = trace), 
        lower = c(0, rep(-Inf, length(inits) - 1)), upper = rep(Inf, 
            length(inits)))
    names(esti$par) <- c("Disp", colnames(X))
    AIC <- 2 * esti$objective + 2 * length(esti$par)
    var <- NULL
    if (est.var) {
        print("Calculating the variance of the estiamtes")
        var <- solve(nH2(pt = esti$par, fun = my.fun))
        colnames(var) <- rownames(var) <- names(esti$par)
    }
    print("Calculating fitted values and residuals")
    PIT <- Qres <- double(nrow(X))
    lam <- exp(X %*% esti$par[-1] + offset)
    lam <- as.double(lam)
    for (ii in 1:nrow(X)) {
        dis <- c(dnbinom(0:(y[ii] - 1), size = esti$par[1], mu = lam[ii], 
            log = FALSE), dnbinom(y[ii], size = esti$par[1], 
            mu = lam[ii], log = FALSE)/2)
        PIT[ii] <- sum(dis[0:y[ii]])
    }
    Qres <- qnorm(PIT)
    print("Exiting")
    out <- list(coef = esti$par, vcov = var, logl = -esti$objective, 
        AIC = AIC, fitted = lam, residuals = Qres, mean.X = mean.X, 
        sd.X = sd.X, formula = NB.form)
    class(out) <- "negBinMod"
    return(out)
}


"nH2" <-
function (pt, fun, accur = c(4, 4), type = "H.Diag", ...) 
{
    H.n <- length(pt)
    derivs <- function(d.x0, ...) {
        nd2(x0 = d.x0, f = fun, m = 1, D.accur = accur[2], ...)
    }
    Hes <- nd2(x0 = pt, f = derivs, D.accur = accur[2], ...)
    Hes <- matrix(Hes, nrow = length(pt))
    Hes <- (Hes + t(Hes))/2
    if (type == "H.Diag") {
        macheps <- .Machine$double.eps
        H.h <- macheps^(1/4) * abs(pt)
        H.f0 <- fun(pt, ...)
        H.m <- length(H.f0)
        if (accur[1] == 2) {
            H.w <- tcrossprod(rep(1, H.m), c(1, -2, 1))
            H.co <- c(-1, 0, 1)
        }
        else {
            H.w <- tcrossprod(rep(1, H.m), c(-1/12, 4/3, -5/2, 
                4/3, -1/12))
            H.co <- c(-2, -1, 0, 1, 2)
        }
        H.n.c <- length(H.co)
        Hes.diag <- double(length = H.n)
        for (ii in 1:H.n) {
            H.temp.f <- matrix(0, H.m, H.n.c)
            for (jj in 1:H.n.c) {
                if (H.co[jj] != 0) {
                  H.xd <- pt + H.h[ii] * H.co[jj] * (1:H.n == 
                    ii)
                  H.temp.f[, jj] <- fun(H.xd, ...)
                }
                else H.temp.f[, jj] <- H.f0
            }
            Hes.diag[ii] <- rowSums(H.w * H.temp.f)/(H.h[ii]^2)
        }
        diag(Hes) <- Hes.diag
    }
    return(Hes)
}


".onLoad" <-
function (libname, pkgname) 
{
    dll.path <- file.path(libname, pkgname, "libs")
    if (nzchar(subarch <- .Platform$r_arch)) 
        dll.path <- file.path(dll.path, subarch)
    this.ext <- paste(sub(".", "[.]", .Platform$dynlib.ext, fixed = TRUE), 
        "$", sep = "")
    dlls <- dir(dll.path, pattern = this.ext, full.names = FALSE)
    names(dlls) <- dlls
    if (length(dlls)) 
        lapply(dlls, function(x) library.dynam(sub(this.ext, 
            "", x), package = pkgname, lib.loc = libname))
}


"predict.MDMnMod" <-
function (object, new.obs, N = NA, S = NA, ...) 
{
    model.nij <- object
    if (S < 2 | is.na(S)) 
        return(list(deriv.eta = NA, nij = NA))
    new.obs <- as.matrix(data.frame(intercept = 1, model.frame(model.nij$formula[-2], 
        new.obs, na.action = NULL)))
    if (!is.na(model.nij$mean.X)[1]) 
        new.obs[-1] <- (new.obs[-1] - model.nij$mean.X)/model.nij$sd.X
    tau <- NA
    cnttau <- 1
    while (is.na(tau) & cnttau < 20) {
        parm.nij <- rmvnorm(1, model.nij$coef, model.nij$vcov)
        tau <- -(new.obs %*% parm.nij[3:length(parm.nij)])
        cnttau <- cnttau + 1
    }
    deriv.j.sampled <- tau[1]
    tau <- tau * log(1:S)
    pij.sampled <- exp(tau)/sum(exp(tau))
    list(deriv.eta = deriv.j.sampled, nij = pij.sampled * N)
}


"predict.negBinMod" <-
function (object, new.obs, offset = 1, ...) 
{
    model.N <- object
    N <- NA
    cntN <- 1
    new.obs <- model.frame(model.N$formula[-2], new.obs, na.action = NULL)
    if (!is.na(model.N$mean.X)[1]) 
        new.obs[-length(new.obs)] <- (new.obs[-length(new.obs)] - 
            model.N$mean.X)/model.N$sd.X
    new.obs <- data.frame(intercept = 1, new.obs[-length(new.obs)])
    while (is.na(N) & cntN < 20) {
        parm.n <- rmvnorm(1, model.N$coef, model.N$vcov)
        tau <- sum(new.obs * parm.n[2:length(parm.n)]) + log(offset)
        N <- rnbinom(1, size = parm.n[1], mu = exp(tau))
        cntN <- cntN + 1
    }
    expect.N <- exp(tau)
    list(N = N, exp.N = expect.N)
}


"predict.truncMod" <-
function (object, new.obs, N = NA, offset = 1, dist = "NB", ...) 
{
    model.S <- object
    if (N < 2 | is.na(N)) 
        return(data.frame(S = NA, expect.S = NA))
    parm.s <- rmvnorm(1, model.S$coef, model.S$vcov)
    new.obs <- model.frame(model.S$formula[-2], new.obs, na.action = NULL)
    if (!is.na(model.S$mean.X)[1]) 
        new.obs[-length(new.obs)] <- (new.obs[-length(new.obs)] - 
            model.S$mean.X)/model.S$sd.X
    new.obs <- data.frame(intercept = 1, new.obs[-length(new.obs)])
    if (dist == "NB") {
        cntS <- 1
        while (parm.s[1] <= 0 & cntS < 20) {
            parm.s <- rmvnorm(1, model.S$coef, model.S$vcov)
            cntS <- cntS + 1
        }
        tau <- sum(new.obs * parm.s[2:length(parm.s)]) + log(offset)
        unscaled.dist <- dnbinom(0:N, size = parm.s[1], mu = exp(tau), 
            log = FALSE)
    }
    if (dist == "Poisson") {
        tau <- sum(new.obs * parm.s) + log(offset)
        unscaled.dist <- dpois(0:N, lambda = exp(tau), log = FALSE)
    }
    scaled.dist <- unscaled.dist/sum(unscaled.dist)
    my.expect <- sum((0:N) * scaled.dist)
    CDF <- c(0, cumsum(scaled.dist))
    U <- runif(1, 0, 1)
    my.which <- tail(which(CDF < U), 1)
    my.S <- (0:N)[my.which]
    list(S = my.S, expect.S = my.expect)
}


"sumLogl" <-
function (tau, theta = NULL, nu = NULL, n, X, logij, ID, dist) 
{
    pij <- exp(-(X %*% tau) * logij)
    LSMN <- double(length = length(unique(ID)))
    kount <- 1
    if (dist == "Mn") {
        for (ii in unique(ID)) {
            pijID <- pij[ID == ii]/sum(pij[ID == ii])
            LSMN[kount] <- dmultinom(x = n[ID == ii], prob = pijID, 
                log = TRUE)
            kount <- kount + 1
        }
        return(sum(LSMN))
    }
    if (dist == "DMn") {
        for (ii in unique(ID)) {
            pijID <- pij[ID == ii]/sum(pij[ID == ii])
            nID <- n[ID == ii]
            LSMN[kount] <- dlDirichletMultinomial(nID, pijID, 
                theta)
            kount <- kount + 1
        }
        return(sum(LSMN))
    }
    if (dist == "MDMn") {
        for (ii in unique(ID)) {
            pijID <- pij[ID == ii]/sum(pij[ID == ii])
            alpha <- theta * pijID
            A <- theta - cumsum(alpha)
            A[A < 0] <- 0
            N <- rev(cumsum(rev(n[ID == ii])))
            S <- length(n[ID == ii])
            LSMN[kount] <- sum(.Call("estBetaBinomial", N, n[ID == 
                ii], alpha, A, 1, -nu, PACKAGE = "RAD")[-S])
            kount <- kount + 1
        }
        return(sum(LSMN))
    }
}


"TLogl" <-
function (tau, od, X, y, offset, t, dist) 
{
    lp <- X %*% tau + offset
    lam <- exp(lp)
    summy <- rep(0, nrow(X))
    if (dist == "NB") {
        lnb <- dnbinom(y, size = od, mu = lam, log = TRUE)
        for (ii in 1:nrow(X)) summy[ii] <- sum(dnbinom(0:t[ii], 
            size = od, mu = lam[ii], log = FALSE))
    }
    if (dist == "Poisson") {
        lnb <- dpois(y, lam, log = TRUE)
        for (ii in 1:nrow(X)) {
            summy[ii] <- sum(dpois(0:t[ii], lam[ii], log = FALSE))
        }
    }
    res <- sum(lnb - log(summy))
    return(res)
}


"TLogl.c" <-
function (tau, od, X, y, offset, t, dist) 
{
    if (dist == "NB") {
        res <- .Call("estllTrunc", as.numeric(y), od, tau, X, 
            offset, t, as.integer(1), PACKAGE = "RAD")
        return(res)
    }
    if (dist == "Poisson") {
        lp <- X %*% tau + offset
        lam <- exp(lp)
        summy <- rep(0, nrow(X))
        lnb <- dpois(y, lam, log = TRUE)
        for (ii in 1:nrow(X)) {
            summy[ii] <- sum(dpois(0:t[ii], lam[ii], log = FALSE))
        }
        res <- sum(lnb - log(summy))
    }
    return(res)
}


"TLogl.old" <-
function (tau, od, X, y, offset, t, dist) 
{
    lp <- X %*% tau + offset
    lam <- exp(lp)
    summy <- rep(0, nrow(X))
    if (dist == "NB") {
        lnb <- dnbinom(y, size = od, mu = lam, log = TRUE)
        for (ii in 1:nrow(X)) summy[ii] <- sum(dnbinom(0:t[ii], 
            size = od, mu = lam[ii], log = FALSE))
    }
    if (dist == "Poisson") {
        lnb <- dpois(y, lam, log = TRUE)
        for (ii in 1:nrow(X)) {
            summy[ii] <- sum(dpois(0:t[ii], lam[ii], log = FALSE))
        }
    }
    res <- sum(lnb - log(summy))
    return(res)
}


"truncMod" <-
function (trunc.form, trunc.pts, data, dist = "NB", scale.covar = FALSE, 
    est.var = TRUE, trace = TRUE) 
{
    temp <- model.frame(trunc.form, as.data.frame(data))
    y <- model.response(temp)
    t <- trunc.pts
    names(y) <- NULL
    X <- model.matrix(trunc.form, data)
    mean.X <- sd.X <- NA
    if (scale.covar & ncol(X) > 1) {
        X.t <- data.frame(X)
        mean.X <- apply(X.t[, 2:ncol(X.t)], 2, mean)
        names(mean.X) <- names(X.t)[2:ncol(X.t)]
        sd.X <- apply(X.t[, 2:ncol(X.t)], 2, sd)
        X[, 2:ncol(X)] <- scale(X[, 2:ncol(X)])
    }
    offset <- model.offset(temp)
    if (is.null(offset)) 
        offset <- rep(0, nrow(temp))
    print("Estimating parameters")
    print("Finding initial values using untruncated model and profile quasi-likelihood")
    if (dist == "Poisson") {
        inits <- glm(trunc.form, data = data, family = poisson(link = "log"))$coef
        my.fun <- function(x) {
            -TLogl(x, 1, X, y, offset, t, dist)
        }
        cat("Ite:     -logl    : ", paste(colnames(X), collapse = " :   "), 
            "\n")
        nam <- colnames(X)
        esti <- nlminb(start = inits, my.fun, control = list(trace = trace))
    }
    if (dist == "NB") {
        fm.nb <- glm.nb(trunc.form, data = data, link = log, 
            control = glm.control(maxit = 50))
        inits <- c(fm.nb$theta, fm.nb$coef)
        my.fun <- function(x) {
            -TLogl(x[-1], x[1], X, y, offset, t, dist)
        }
        cat("Ite:     -logl    :  Disp    ", paste(colnames(X), 
            collapse = " :   "), "\n")
        nam <- c("theta", colnames(X))
        esti <- nlminb(start = inits, my.fun, control = list(trace = trace), 
            lower = c(0, rep(-Inf, length(inits) - 1)), upper = rep(Inf, 
                length(inits)))
    }
    names(esti$par) <- nam
    AIC <- 2 * esti$objective + 2 * length(esti$par)
    var <- 0
    if (est.var) {
        print("Calculating the variance of the estiamtes")
        var <- solve(nH2(pt = esti$par, fun = my.fun))
        colnames(var) <- rownames(var) <- names(esti$par)
    }
    print("Calculating Fitted Values and Residuals")
    expect <- expect2 <- sds <- pearson <- sums <- PIT <- Qres <- double(nrow(X))
    if (dist == "Poisson") {
        lam <- as.double(exp(X %*% esti$par + offset))
        for (ii in 1:nrow(X)) {
            dis <- dpois(0:t[ii], lam[ii], log = FALSE)
            sums[ii] <- sum(dis)
            dis <- dis/sums[ii]
            PIT[ii] <- sum(dis[0:(y[ii] - 1)]) + dis[y[ii]]/2
            expect[ii] <- sum(0:t[ii] * dis)
            expect2[ii] <- sum((0:t[ii])^2 * dis)
        }
    }
    if (dist == "NB") {
        lam <- as.double(exp(X %*% esti$par[-1] + offset))
        for (ii in 1:nrow(X)) {
            dis <- dnbinom(0:t[ii], size = esti$par[1], mu = lam[ii], 
                log = FALSE)
            sums[ii] <- sum(dis)
            dis <- dis/sums[ii]
            PIT[ii] <- sum(dis[0:(y[ii] - 1)]) + dis[y[ii]]/2
            expect[ii] <- sum(0:t[ii] * dis)
            expect2[ii] <- sum((0:t[ii])^2 * dis)
        }
    }
    sds <- sqrt(expect2 - (expect)^2)
    Qres <- qnorm(PIT)
    print("Exiting")
    out <- list(coef = esti$par, vcov = var, logl = -esti$objective, 
        AIC = AIC, residuals = Qres, fitted = expect, sds = sds, 
        sums = sums, mean.X = mean.X, sd.X = sd.X, formula = trunc.form)
    class(out) <- "truncMod"
    return(out)
}

