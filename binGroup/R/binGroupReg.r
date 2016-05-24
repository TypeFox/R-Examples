gtreg <- function(formula, data, groupn, retest = NULL, sens = 1, spec = 1,
        linkf = c("logit", "probit", "cloglog"),
        method = c("Vansteelandt", "Xie"), sens.ind = NULL, spec.ind = NULL,
        start = NULL, control = gt.control(...), ...) {

    call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "groupn"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    gr <- model.extract(mf, "groupn")
    if (!is.na(pos <- match(deparse(substitute(retest)), names(data))))
        retest <- data[, pos]
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm))
            names(Y) <- nm
    }
    X <- if (!is.empty.model(mt))
        model.matrix(mt, mf)
    else matrix(, NROW(Y), 0)
    linkf <- match.arg(linkf)
    if ((method <- match.arg(method)) == "Vansteelandt") {
        if (!is.null(retest))
            warning("Retests cannot be used with Vansteelandt's method.")
        fit <- gtreg.fit(Y, X, gr, sens, spec, linkf, start)
    }
    else {
        if (is.null(retest)) 
            fit <- EM(Y, X, gr, sens, spec, linkf, start, control)
        else fit <-  EM.ret(Y, X, gr, retest, sens, spec, linkf,
             sens.ind, spec.ind, start, control)
    }
    fit <- c(fit, list(call = call, formula = formula, method = method,
          link = linkf, terms = mt))
    class(fit) <- "gt"
    fit

}


gtreg.fit <- function (Y, X, groupn, sens, spec, linkf, start = NULL)
{
    z <- tapply(Y, groupn, tail, n = 1)
    num.g <- max(groupn)
    K <- ncol(X)
    sam <- length(Y)
    if (is.null(start)) {
        if (K == 1) {
            cova.mean <- as.matrix(tapply(X, groupn, mean)) 
            optim.meth <- "BFGS"
        } 
        else {
            temp <- by(X, groupn, colMeans)
            cova.mean <- do.call(rbind, temp)
            optim.meth <- "Nelder-Mead"
        }
        beta.group <- glm.fit(cova.mean, as.vector(z), 
            family = binomial(link = linkf))$coefficients
    } 
    else {
        beta.group <- start
        names(beta.group) <- dimnames(X)[[2]]
        optim.meth <- ifelse(K == 1, "BFGS", "Nelder-Mead")
    }   
    logL <- function(beta) {
        pijk <- switch(linkf, logit = plogis(X %*% beta), probit = pnorm(X %*%
            beta), cloglog = 1 - exp(-exp(X %*% beta)))
        prodp <- tapply(1 - pijk, groupn, prod)
        -sum(z * log(sens + (1 - sens - spec) * prodp) +
            (1 - z) * log(1 - sens - (1 - sens - spec) * prodp))
    }
    mod.fit <- optim(par = beta.group, fn = logL, method = optim.meth,
        control = list(trace = 0, maxit = 1000), hessian = TRUE)
    if (det(mod.fit$hessian) == 0)
        mod.fit <- optim(par = beta.group, fn = logL, method = "SANN", hessian = TRUE)
    logL0 <- function(beta) {
        inter <- rep(beta, sam)
        pijk <- switch(linkf, logit = plogis(inter), probit = pnorm(inter),
            cloglog = 1 - exp(-exp(inter)))
        prodp <- tapply(1 - pijk, groupn, prod)
        -sum(z * log(sens + (1 - sens - spec) * prodp) +
            (1 - z) * log(1 - sens - (1 - sens - spec) * prodp))
    }
    mod.fit0 <- optim(par = binomial()$linkfun(mean(z)), 
        fn = logL0, method = "BFGS", control = list(trace = 0, maxit = 1000))
    nulld <- 2 * mod.fit0$value
    residd <- 2 * mod.fit$value
    xib <- X %*% mod.fit$par
    pijk <- switch(linkf, logit = plogis(xib), probit = pnorm(xib),
        cloglog = 1 - exp(-exp(xib)))
    prodp <- tapply(1 - pijk, groupn, prod)
    zhat <- sens + (1 - sens - spec) * prodp
    residual <- z - zhat
    aic <- residd + 2 * K
    if (mod.fit$convergence == 0)
        counts <- mod.fit$counts[[1]]
    else warning("Maximum number of iterations exceeded.")
    list(coefficients = mod.fit$par, hessian = mod.fit$hessian,
        fitted.values = zhat, deviance = residd, df.residual = num.g - K,
        null.deviance = nulld, df.null = num.g - 1, aic = aic, counts = counts,
        residuals = residual, z = z)
}


EM <- function (Y, X, groupn, sens, spec, linkf, start = NULL, control = gt.control())
{
    if (control$time)
        start.time <- proc.time()
    z <- tapply(Y, groupn, tail, n = 1)
    num.g <- max(groupn)
    K <- ncol(X)
    if (is.null(start)) {
        if (K == 1)
            cova.mean <- as.matrix(tapply(X, groupn, mean))
        else {
            temp <- by(X, groupn, colMeans)
            cova.mean <- do.call(rbind, temp)
        }
        beta.old <- lm.fit(cova.mean, z)$coefficients
    } 
    else beta.old <- start
    sam <- length(Y)
    vec <- 1:sam
    group.sizes <- tapply(Y, groupn, length)
    diff <- 1
    counts <- 1
    extra.loop <- FALSE
    next.loop <- TRUE
    while (next.loop) {
        xib <- X %*% beta.old
        pijk <- switch(linkf, logit = plogis(xib),
            probit = pnorm(xib), cloglog = 1 - exp(-exp(xib)))
        prodp <- tapply(1 - pijk, groupn, prod)
        den <- rep((1 - spec) * prodp + sens * (1 - prodp), group.sizes)
        den2 <- rep(spec * prodp + (1 - sens) * (1 - prodp),
            group.sizes)
        expect <- rep(NA, times = sam)
        for (i in vec) {
            if (Y[i] == 0)
                expect[i] <- (1 - sens) * pijk[i]/den2[i]
            else expect[i] <- sens * pijk[i]/den[i]
        }
        if (!extra.loop) {
            suppress <- function(w) 
                if(any(grepl("non-integer #successes in a binomial glm", w))) 
                   invokeRestart("muffleWarning")
            mod.fit <- withCallingHandlers(glm.fit(X, expect, 
                family = binomial(link = linkf)), warning = suppress)
            diff <- max(abs((beta.old - mod.fit$coefficients)/beta.old))
            beta.old <- mod.fit$coefficients
            if (control$trace)
                cat("beta is", beta.old, "\tdiff is", diff, "\n")
            counts <- counts + 1
            if (diff <= control$tol || counts > control$maxit) 
                extra.loop <- TRUE
        } 
        else next.loop <- FALSE
    }    
    erf <- 2 * pijk - 1
    pt1 <- switch(linkf, logit = -exp(xib)/(1 + exp(xib))^2,
        probit = sqrt(2) * xib * exp(-xib^2/2)/(sqrt(pi) * (1 -
            erf)) - 2 * exp(-xib^2)/(pi * (1 - erf)^2), cloglog = -exp(xib))
    pt2 <- switch(linkf, logit = 0, probit = (8 * exp(-xib^2/2) *
        erf + 2 * xib * sqrt(2 * pi) * erf^2 - 2 * xib * sqrt(2 *
        pi)) * exp(-xib^2/2)/((1 + erf)^2 * pi * (1 - erf)^2),
        cloglog = -(exp(xib - exp(xib)) + exp(2 * xib - exp(xib)) -
            exp(xib))/(exp(-exp(xib)) - 1)^2)
    nm <- pt1 + expect * pt2
    sign1 <- as.vector(sign(nm))
    nn <- as.vector(sqrt(abs(nm)))
    x2 <- X * nn
    m <- (t(x2) %*% (sign1 * x2))
    b <- array(NA, c(K, K, sum(group.sizes^2)))
    p <- 1
    for (i in vec) for (j in vec[groupn == groupn[i]]) {
        wii <- ifelse(i == j, expect[i] - expect[i]^2, expect[i] *
            (pijk[j] - expect[j]))
        coe <- switch(linkf, logit = 1, probit = 8 * exp(-(xib[i]^2 +
            xib[j]^2)/2)/((1 - erf[i]^2) * (1 - erf[j]^2) * pi),
            cloglog = exp(xib[i] + xib[j])/((exp(-exp(xib[i])) -
                1) * (exp(-exp(xib[j])) - 1)))
        b[, , p] <- wii * coe * X[i, ] %*% t(X[j, ])
        p <- p + 1
    }
    m1 <- apply(b, c(1, 2), sum)
    H <- -(m + m1)
    zhat <- sens + (1 - sens - spec) * prodp
    residual <- z - zhat
    residd <- -2 * sum(z * log(zhat) + (1 - z) * log(1 - zhat))
    logL0 <- function(beta) {
        inter <- rep(beta, sam)
        pijk <- switch(linkf, logit = plogis(inter), probit = pnorm(inter),
            cloglog = 1 - exp(-exp(inter)))
        prodp <- tapply(1 - pijk, groupn, prod)
        -sum(z * log(sens + (1 - sens - spec) * prodp) +
            (1 - z) * log(1 - sens - (1 - sens - spec) * prodp))
    }
    mod.fit0 <- optim(par = binomial()$linkfun(mean(z)), fn = logL0,
        method = "BFGS", control = list(trace = 0, maxit = 1000))
    nulld <- 2 * mod.fit0$value
    aic <- residd + 2 * K
    if (diff > control$tol && counts > control$maxit)
        warning("EM algorithm did not converge.")
    if (control$time) {
        end.time <- proc.time()
        save.time <- end.time - start.time
        cat("\n Number of minutes running:", round(save.time[3]/60, 2), "\n \n")
    }
    list(coefficients = beta.old, hessian = H, fitted.values = zhat,
        deviance = residd, df.residual = num.g - K, null.deviance = nulld,
        df.null = num.g - 1, aic = aic, counts = counts - 1, residuals = residual,
        z = z)
}


EM.ret <- function (Y, X, groupn, ret, sens, spec, linkf,
    sens.ind, spec.ind,
    start = NULL, control = gt.control())
{
    if (control$time)
        start.time <- proc.time()
    if (is.null(sens.ind))
        sens.ind <- sens
    if (is.null(spec.ind))
        spec.ind <- spec
    z <- tapply(Y, groupn, tail, n = 1)
    num.g <- max(groupn)
    K <- ncol(X)
    if (is.null(start)) {
        if (K == 1)
            cova.mean <- as.matrix(tapply(X, groupn, mean))
        else {
            temp <- by(X, groupn, colMeans)
            cova.mean <- do.call(rbind, temp)
        }
        beta.old <- lm.fit(cova.mean, z)$coefficients
    }
    else beta.old <- start
    sam <- length(Y)
    vec <- 1:sam
    group.sizes <- tapply(Y, groupn, length)
    diff <- 1
    counts <- 1
    extra.loop <- FALSE
    next.loop <- TRUE
    a0 <- ifelse(ret == 1, sens.ind, 1 - sens.ind)
    a1 <- ifelse(ret == 0, spec.ind, 1 - spec.ind)
    while (next.loop) {
        xib <- X %*% beta.old
        pijk <- switch(linkf, logit = plogis(xib),
            probit = pnorm(xib), cloglog = 1 -
                exp(-exp(xib)))
        erf <- 2 * pijk - 1
        prodp <- tapply(1 - pijk, groupn, prod)
        den2 <- rep(spec * prodp + (1 - sens) * (1 - prodp),
            group.sizes)
        expect <- rep(NA, times = sam)
        i <- 1
        while (i <= sam) {
            if (Y[i] == 0)
                expect[i] <- (1 - sens) * pijk[i]/den2[i]
            else {
                vec1 <- vec[groupn == groupn[i]]
                mb2 <- 1
                for (l in vec1) {
                     temp <- a0[l] * pijk[l] + a1[l] * (1 - pijk[l])
                     mb2 <- mb2 * temp
                }
                null <- 1
                for (l in vec1) {
                     temp <- a1[l] * (1 - pijk[l])
                     null <- null * temp
                }                                                       
                den <- mb2 * sens + null * (1 - sens - spec)
                for (l1 in vec1) {                         
                     temp <- a0[l1] * pijk[l1] + a1[l1] * (1 - pijk[l1])
                     num <- mb2/temp * a0[l1] * pijk[l1] * sens
                     expect[l1] <- num/den
                }
                i <- l1
            }
            i <- i + 1
        }
        expect[expect > 1] <- 1
        expect[expect < 0] <- 0
        if (!extra.loop) {
            suppress <- function(w)
                if (any(grepl("non-integer #successes in a binomial glm", w)))
                   invokeRestart("muffleWarning")
            mod.fit <- withCallingHandlers(glm.fit(X, expect,
                family = binomial(link = linkf)), warning = suppress)
            diff <- max(abs((beta.old - mod.fit$coefficients)/beta.old))
            beta.old <- mod.fit$coefficients
            if (control$trace)
                cat("beta is", beta.old, "\tdiff is", diff, "\n")
            counts <- counts + 1
            if (diff <= control$tol || counts > control$maxit) 
                extra.loop <- TRUE
        } 
        else next.loop <- FALSE
    }
    pt1 <- switch(linkf, logit = -exp(xib)/(1 + exp(xib))^2,
        probit = sqrt(2) * xib * exp(-xib^2/2)/(sqrt(pi) * (1 -
            erf)) - 2 * exp(-xib^2)/(pi * (1 - erf)^2), cloglog = -exp(xib))
    pt2 <- switch(linkf, logit = 0, probit = (8 * exp(-xib^2/2) *
        erf + 2 * xib * sqrt(2 * pi) * erf^2 - 2 * xib * sqrt(2 *
        pi)) * exp(-xib^2/2)/((1 + erf)^2 * pi * (1 - erf)^2),
        cloglog = -(exp(xib - exp(xib)) + exp(2 * xib - exp(xib)) -
            exp(xib))/(exp(-exp(xib)) - 1)^2)
    nm <- pt1 + expect * pt2
    sign1 <- as.vector(sign(nm))
    nn <- as.vector(sqrt(abs(nm)))
    x2 <- X * nn
    m <- (t(x2) %*% (sign1 * x2))
    m1 <- 0    
    for (i in vec) {
        vec1 <- vec[groupn == groupn[i]]
        if (Y[i] == 0) {
            for (j in vec1) {
                 coe <- switch(linkf, logit = 1, probit = 8 * exp(-(xib[i]^2 +
                     xib[j]^2)/2)/((1 - erf[i]^2) * (1 - erf[j]^2) * pi),
                     cloglog = exp(xib[i] + xib[j])/((exp(-exp(xib[i])) -
                     1) * (exp(-exp(xib[j])) - 1)))            
                 wii <- ifelse(i == j, expect[i] - expect[i]^2, expect[i] *
                    (pijk[j] - expect[j]))
                 tim <- wii * coe * X[i, ] %*% t(X[j, ])
                 m1 <- m1 + tim
            }                
        }        
        else {         
            for (j in vec1) {      
                 temp <- a0[j] * pijk[j] + a1[j] * (1 - pijk[j])
                 eii <- expect[i]/temp * a0[j] * pijk[j]
                 wii <- ifelse(i == j, expect[i] - expect[i]^2, eii - expect[i] * expect[j])
                 coe <- switch(linkf, logit = 1, probit = 8 * exp(-(xib[i]^2 +
                      xib[j]^2)/2)/((1 - erf[i]^2) * (1 - erf[j]^2) * pi),
                      cloglog = exp(xib[i] + xib[j])/((exp(-exp(xib[i])) -
                        1) * (exp(-exp(xib[j])) - 1)))
                 tim <- wii * coe * X[i, ] %*% t(X[j, ])
                 m1 <- m1 + tim
            }
        }
    }
    H <- -(m + m1)
    zhat <- sens + (1 - sens - spec) * prodp
    residual <- z - zhat
    logl <- 0
    for (grn in 1:num.g) {
         if (z[grn] == 1) {
             vec1 <- vec[groupn == grn]
             mb2 <- 1
             for (l in vec1) {
                  temp <- a0[l] * pijk[l] + a1[l] * (1 - pijk[l])
                  mb2 <- mb2 * temp
             }
             null <- 1
             for (l in vec1) {
                  temp <- a1[l] * (1 - pijk[l])
                  null <- null * temp
             }                                                       
             prob1 <- mb2 * sens + null * (1 - sens - spec)
         } else prob1 <- 1 - zhat[grn]
         logl <- logl - log(prob1)
    }
    aic <- 2 * logl + 2 * K
    if (diff > control$tol && counts > control$maxit)
        warning("EM algorithm did not converge.")
    if (control$time) {
        end.time <- proc.time()
        save.time <- end.time - start.time
        cat("\n Number of minutes running:", round(save.time[3]/60, 2), "\n \n")
    }
    list(coefficients = beta.old, hessian = H, fitted.values = zhat,
        deviance = 2 * logl, aic = aic, counts = counts - 1, residuals = residual, 
        z = z)
}


gt.control <- function (tol = 0.0001, n.gibbs = 1000, n.burnin = 20, 
    maxit = 500, trace = FALSE, time = TRUE) 
{
    if (!is.numeric(tol) || tol <= 0) 
        stop("value of 'tol' must be > 0")
    if (round(n.gibbs) != n.gibbs || n.gibbs <= 0) 
        stop("value of 'n.gibbs' must be a positive integer")
    if (round(n.burnin) != n.burnin || n.burnin <= 0) 
        stop("value of 'n.burnin' must be a positive integer")    
    if (!is.numeric(maxit) || maxit <= 0) 
        stop("maximum number of iterations must be > 0")
    list(tol = tol, n.gibbs = n.gibbs, n.burnin = n.burnin, maxit = maxit, 
        trace = trace, time = time)
}


gtreg.mp <- function (formula, data, coln, rown, arrayn, retest = NULL,
    sens = 1, spec = 1, linkf = c("logit", "probit", "cloglog"), 
    sens.ind = NULL, spec.ind = NULL, start = NULL, control = gt.control(...), ...)
{
    call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "coln", "rown",
        "arrayn"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    arrayn <- model.extract(mf, "arrayn")
    rown <- model.extract(mf, "rown")
    coln <- model.extract(mf, "coln")
    if (!is.na(pos <- match(deparse(substitute(retest)), names(data))))
        retest <- data[, pos]
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm))
            names(Y) <- nm
    }
    X <- if (!is.empty.model(mt))
        model.matrix(mt, mf)
    else matrix(, NROW(Y), 0)
    linkf <- match.arg(linkf)
    fit <- EM.mp(Y[, 1], Y[, 2], X, coln, rown, arrayn, retest,
        sens, spec, linkf, sens.ind, spec.ind, start, control)
    fit <- c(fit, list(call = call, formula = formula, link = linkf,
        terms = mt))
    class(fit) <- c("gt.mp", "gt")
    fit
}


EM.mp <- function (col.resp, row.resp, X, coln, rown, sqn, ret, sens,
    spec, linkf, sens.ind, spec.ind, start = NULL, control = gt.control())
{
    if (control$time)
        start.time <- proc.time()
    if (is.null(sens.ind))
        sens.ind <- sens
    if (is.null(spec.ind))
        spec.ind <- spec
    len <- max(sqn)
    diff <- 1
    counts <- 1
    sam <- length(sqn)
    col.groupn <- coln[sqn == 1]
    if (len > 1) {
        for (i in 2:len) {
            temp <- max(col.groupn) + coln[sqn == i]
            col.groupn <- c(col.groupn, temp)
        }
    }
    if (is.null(start)) {
        mod.fit <- try(gtreg.fit(col.resp, X, col.groupn,
            sens, spec, linkf))
        if (class(mod.fit) == "try-error") {
            row.groupn <- rown[sqn == 1]
            if (len > 1) {
                for (i in 2:len) {
                  temp <- max(row.groupn) + rown[sqn == i]
                  row.groupn <- c(row.groupn, temp)
                }
            }
            mod.fit <- gtreg.fit(row.resp, X, row.groupn,
                sens, spec, linkf)
        }
        beta.old <- mod.fit$coefficients
    }
    else beta.old <- start
    extra.loop <- FALSE
    next.loop <- TRUE
    while (next.loop) {
        xib <- X %*% beta.old
        pijk.all <- switch(linkf, logit = plogis(xib),
            probit = pnorm(xib), cloglog = 1 - exp(-exp(xib)))
        expect.all <- numeric(0)
        mat2 <- index <- 0        
        erf <- 2 * pijk.all - 1
        for (arrayn in 1:len) {
            index.r <- index.c <- vector("logical", length = sam)
            for (i in 1:sam) {
                if (rown[i] == 1 && sqn[i] == arrayn)
                    index.c[i] <- TRUE
                else index.c[i] <- FALSE
                if (coln[i] == 1 && sqn[i] == arrayn)
                    index.r[i] <- TRUE
                else index.r[i] <- FALSE
            }
            n.row <- max(rown[index.r])
            n.col <- max(coln[index.c])
            rowresp <- row.resp[index.r]
            colresp <- col.resp[index.c]
            index <- max(index) + 1:(n.row * n.col)
            if (!is.null(ret)) {
                re.ind <- na.omit(cbind(coln[sqn == arrayn],
                  rown[sqn == arrayn], ret[sqn == arrayn]))
                re <- ifelse(re.ind[, 3] == 1, sens.ind, 1 -
                  sens.ind)
                re1 <- ifelse(re.ind[, 3] == 0, spec.ind, 1 -
                  spec.ind)
            }
            pijk <- matrix(pijk.all[sqn == arrayn], nrow = n.row)
            a <- ifelse(rowresp == 1, sens, 1 - sens)
            b <- ifelse(colresp == 1, sens, 1 - sens)
            a1 <- ifelse(rowresp == 0, spec, 1 - spec)
            b1 <- ifelse(colresp == 0, spec, 1 - spec)
            mat <- array(NA, c(n.row, n.col, control$n.gibbs))
            y <- matrix(0, nrow = n.row, ncol = n.col)
            for (k in 1:(control$n.gibbs + control$n.burnin)) {
                l <- 1
                for (j in 1:n.col) for (i in 1:n.row) {
                  num <- a[i] * b[j] * pijk[i, j]
                  den.r <- ifelse(sum(y[i, ]) - y[i, j] > 0,
                    a[i], a1[i])
                  den.c <- ifelse(sum(y[, j]) - y[i, j] > 0,
                    b[j], b1[j])
                  den2 <- den.r * den.c * (1 - pijk[i, j])
                  if (!is.null(ret)) {
                    if (l <= length(re) && j == re.ind[l, 1] &&
                      i == re.ind[l, 2]) {
                      num <- num * re[l]
                      den2 <- den2 * re1[l]
                      l <- l + 1
                    }
                  }
                  den <- num + den2
                  if (den != 0) {
                    cond.p <- num/den
                    y[i, j] <- rbinom(1, 1, cond.p)
                  }
                  else y[i, j] <- 0
                }
                if (k > control$n.burnin) {
                    mat[, , k - control$n.burnin] <- y
                    vec <- as.vector(y)
                    if (extra.loop)
                        for (i1 in index[vec == 1]) for (j1 in index[vec ==
                          1]) {
                          bq <- switch(linkf, logit = 1, probit = 8 *
                            exp(-(xib[i1]^2 + xib[j1]^2)/2)/((1 - erf[i1]^2) *
                            (1 - erf[j1]^2) * pi), cloglog = exp(xib[i1] +
                            xib[j1])/((exp(-exp(xib[i1])) - 1) * (exp(-exp(xib[j1])) -
                            1))) * X[i1, ] %*% t(X[j1, ])
                          mat2 <- mat2 + bq
                        }
                }
            }
            expect.m <- apply(mat, c(1, 2), mean)
            expect <- as.vector(expect.m)
            expect.all <- c(expect.all, expect)
        }
        if (!extra.loop) {
            suppress <- function(w) 
                if(any(grepl("non-integer #successes in a binomial glm", w))) 
                   invokeRestart("muffleWarning")
            mod.fit <- withCallingHandlers(glm.fit(X, expect.all, 
                family = binomial(link = linkf)), warning = suppress)
            diff <- max(abs((beta.old - mod.fit$coefficients)/beta.old))
            beta.old <- mod.fit$coefficients
            if (control$trace)
                cat("beta is", beta.old, "\tdiff is", diff, "\n")
            counts <- counts + 1
            if (diff <= control$tol || counts > control$maxit) 
                extra.loop <- TRUE
        } 
        else next.loop <- FALSE
    }
    index <- 0
    first <- mat2/control$n.gibbs
    second <- 0
    for (arrayn in 1:len) {
        n.row <- max(rown[sqn == arrayn])
        n.col <- max(coln[sqn == arrayn])
        index <- max(index) + 1:(n.row * n.col)
        expect <- expect.all[index]
        for (i1 in index) for (j1 in index) {
            coe <- switch(linkf, logit = 1, probit = 8 * exp(-(xib[i1]^2 +
                xib[j1]^2)/2)/((1 - erf[i1]^2) * (1 - erf[j1]^2) *
                pi), cloglog = exp(xib[i1] + xib[j1])/((exp(-exp(xib[i1])) -
                1) * (exp(-exp(xib[j1])) - 1)))
            tim <- expect.all[i1] * expect.all[j1] * coe * X[i1,
                ] %*% t(X[j1, ])
            second <- second + tim
        }
    }
    m1 <- first - second
    pt1 <- switch(linkf, logit = -exp(xib)/(1 + exp(xib))^2,
        probit = sqrt(2) * xib * exp(-xib^2/2)/(sqrt(pi) * (1 -
            erf)) - 2 * exp(-xib^2)/(pi * (1 - erf)^2), cloglog = -exp(xib))
    pt2 <- switch(linkf, logit = 0, probit = (8 * exp(-xib^2/2) *
        erf + 2 * xib * sqrt(2 * pi) * erf^2 - 2 * xib * sqrt(2 *
        pi)) * exp(-xib^2/2)/((1 + erf)^2 * pi * (1 - erf)^2),
        cloglog = -(exp(xib - exp(xib)) + exp(2 * xib - exp(xib)) -
            exp(xib))/(exp(-exp(xib)) - 1)^2)
    nm <- pt1 + expect.all * pt2
    sign1 <- as.vector(sign(nm))
    nn <- as.vector(sqrt(abs(nm)))
    x2 <- X * nn
    m <- (t(x2) %*% (sign1 * x2))
    H <- -(m + m1)
    if (diff > control$tol && counts > control$maxit)
        warning("EM algorithm did not converge.")
    if (control$time) {
        end.time <- proc.time()
        save.time <- end.time - start.time
        cat("\n Number of minutes running:", round(save.time[3]/60, 2), "\n \n")
    }
    list(coefficients = beta.old, hessian = H, Gibbs.sample.size = control$n.gibbs,
        counts = counts - 1)
}


sim.gt <- function (x = NULL, gshape = 20, gscale = 2, par,
    linkf = c("logit", "probit", "cloglog"),
    sample.size, group.size, sens = 1, spec = 1, sens.ind = NULL, spec.ind = NULL)
{
    if (is.null(sens.ind))
        sens.ind <- sens
    if (is.null(spec.ind))
        spec.ind <- spec
    if (is.null(x)) {
        x <- rgamma(n = sample.size, shape = gshape, scale = gscale)
        X <- cbind(1, x)
    }
    else {
        X <- cbind(1, x)
        sample.size <- nrow(X)
    }
    linkf <- match.arg(linkf)
    pijk <- switch(linkf, logit = plogis(X %*% par),
            probit = pnorm(X %*% par),
            cloglog = 1 - exp(-exp(X %*% par)))
    ind <- rbinom(n = sample.size, size = 1, prob = pijk)
    num.g <- ceiling(sample.size/group.size)
    vec <- 1:sample.size
    groupn <- rep(1:num.g, each = group.size)[vec]
    save.sum <- tapply(ind, groupn, sum)
    save.group <- as.vector(ifelse(save.sum > 0, 1, 0))
    save.obs <- rep(NA, num.g)
    ret <- rep(NA, sample.size)
    for (i in 1:num.g)
         save.obs[i] <- ifelse(save.group[i] == 1, rbinom(1, 1, sens),
                1 - rbinom(1, 1, spec))
    gres <- rep(save.obs, each = group.size)[vec]
    for (i in vec) {
        if (gres[i] == 1)
            ret[i] <- ifelse(ind[i] == 1, rbinom(1,
                 1, sens.ind), 1 - rbinom(1, 1, spec.ind))
    }   
    grd <- data.frame(gres = gres, x = x, groupn = groupn, ind = ind, retest = ret)
    if (ncol(X) > 2)
        for (i in 2:ncol(X))
             colnames(grd)[i] <- paste("x", i - 1, sep="")
    grd
}


sim.mp <- function (x = NULL, gshape = 20, gscale = 2, par,
    linkf = c("logit", "probit", "cloglog"),
    n.row, n.col, sens = 1, spec = 1, sens.ind = NULL, spec.ind = NULL)
{
    if (is.null(sens.ind))
        sens.ind <- sens
    if (is.null(spec.ind))
        spec.ind <- spec
    if (length(n.row) != length(n.col))
        stop("vector n.row and n.col must have the same length")
    linkf <- match.arg(linkf)
    if (is.null(x)) {
        sample.size <- sum(n.col * n.row)
        x <- rgamma(n = sample.size, shape = gshape, scale = gscale)
        X <- cbind(1, x)
    }
    else {
        X <- cbind(1, x)
        sample.size <- nrow(X)
        if (sum(n.col * n.row) != sample.size)
            stop("n.row and n.col not consistent with the sample size")
    }
    len <- length(n.row)
    pijk <- switch(linkf, logit = plogis(X %*% par),
            probit = pnorm(X %*% par),
            cloglog = 1 - exp(-exp(X %*% par)))
    ind <- rbinom(n = sample.size, size = 1, prob = pijk)
    individual <- col.groupn <- row.groupn <- numeric(0)
    rowr <- colr <- numeric(0)
    ret <- rep(NA, sample.size)
    for (i in 1:len) {
        if (i > 1)
            index <- seq(max(index) + 1, length = (n.row * n.col)[i])
        else index <- 1:(n.row * n.col)[1]
        indm <- matrix(ind[index], nrow = n.row[i])
        col.resp <- apply(indm, MARGIN = 2, FUN = sum)
        col.resp <- ifelse(col.resp > 0, 1, 0)
        col.err <- rep(NA, n.col[i])
        for (j in 1:n.col[i])
             col.err[j] <- ifelse(col.resp[j] == 1, rbinom(1, 1, sens),
                                 1 - rbinom(1, 1, spec))
        row.resp <- apply(indm, MARGIN = 1, FUN = sum)
        row.resp <- ifelse(row.resp > 0, 1, 0)
        row.err <- rep(NA, n.row[i])
        for (j in 1:n.row[i])
             row.err[j] <- ifelse(row.resp[j] == 1, rbinom(1, 1, sens),
                                 1 - rbinom(1, 1, spec))
        temp.c <- rep(1:n.col[i], each = n.row[i])
        col.groupn <- c(col.groupn, temp.c)
        temp.r <- rep(1:n.row[i], n.col[i])
        row.groupn <- c(row.groupn, temp.r)
        temp2.c <- rep(col.err, each = n.row[i])
        colr <- c(colr, temp2.c)
        temp2.r <- rep(row.err, n.col[i])
        rowr <- c(rowr, temp2.r)
        if (all(row.err == 0)) {
            for (j in index) {
                 if (colr[j] == 1)
                     ret[j] <- ifelse(ind[j] == 1, rbinom(1,
                        1, sens.ind), 1 - rbinom(1, 1, spec.ind))
            }
        }
        else {
            if (all(col.err == 0)) {
                for (j in index) {
                     if (rowr[j] == 1)
                         ret[j] <- ifelse(ind[j] == 1, rbinom(1,
                            1, sens.ind), 1 - rbinom(1, 1, spec.ind))
                }
            }
            else {
                for (j in index) {
                     if (rowr[j] == 1 && colr[j] == 1)
                         ret[j] <- ifelse(ind[j] == 1, rbinom(1, 1,
                            sens.ind), 1 - rbinom(1, 1, spec.ind))
                }
            }
        }
        individual <- c(individual, list(indm))
    }
    sq <- rep(1:len, n.col * n.row)
    if (all(colr == 0) && all(rowr == 0))
        return(NULL)
    grd <- data.frame(x = x, col.resp = colr,
        row.resp = rowr, coln = col.groupn, rown = row.groupn,
        arrayn = sq, retest = ret)
    if (ncol(X) > 2)
        for (i in 1:(ncol(X) - 1))
             colnames(grd)[i] <- paste("x", i, sep="")
    list(dframe = grd, ind = individual, prob = as.vector(pijk))
}


summary.gt <- function (object, ...)
{
    coef.p <- object$coefficients
    cov.mat <- solve(object$hessian)
    dimnames(cov.mat) <- list(names(coef.p), names(coef.p))
    var.cf <- diag(cov.mat)
    s.err <- sqrt(var.cf)
    zvalue <- coef.p/s.err
    dn <- c("Estimate", "Std. Error")
    pvalue <- 2 * pnorm(-abs(zvalue))
    coef.table <- cbind(coef.p, s.err, zvalue, pvalue)
    dimnames(coef.table) <- list(names(coef.p), c(dn, "z value",
        "Pr(>|z|)"))
    keep <- match(c("call", "link", "aic", "deviance", "df.residual",
        "null.deviance", "df.null", "counts", "method", "z"),
        names(object), 0)
    ans <- c(object[keep], list(coefficients = coef.table, deviance.resid = residuals(object,
        type = "deviance"), cov.mat = cov.mat))
    class(ans) <- "summary.gt"
    return(ans)
}


print.summary.gt <- function (x, digits = max(3, getOption("digits") - 3),
    signif.stars = getOption("show.signif.stars"), ...)
{
    obj <- x
    cat("\nCall:\n")
    cat(paste(deparse(obj$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
    cat("Deviance Residuals: \n")
    if (length(obj$z) > 5) {
        obj$deviance.resid <- quantile(obj$deviance.resid, na.rm = TRUE)
        names(obj$deviance.resid) <- c("Min", "1Q", "Median",
            "3Q", "Max")
    }
    print.default(obj$deviance.resid, digits = digits, na.print = "",
        print.gap = 2)
    cat("\nCoefficients:\n")
    coefs <- obj$coefficients
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
        na.print = "NA", ...)
    if (!is.null(unlist(obj["df.null"])))
        cat("\n", apply(cbind(paste(format(c("Null", "Residual"),
            justify = "right"), "deviance:"), format(unlist(obj[c("null.deviance",
            "deviance")]), digits = 4), " on", format(unlist(obj[c("df.null",
            "df.residual")])), " degrees of freedom\n"), 1, paste,
            collapse = " "), sep = "")
    if (obj$method == "Vansteelandt")
        cat("AIC: ", format(obj$aic, digits = 4), "\n\n", "Number of iterations in optim(): ",
            obj$counts, "\n", sep = "")
    else {
        cat("AIC: ", format(obj$aic, digits = 4), "\n\n", "Number of iterations in EM: ",
            obj$counts, "\n", sep = "")
    }
    cat("\n")
    invisible(obj)
}


predict.gt <- function (object, newdata, type = c("link", "response"), se.fit = FALSE, 
    conf.level = NULL, na.action = na.pass, ...) 
{
    tt <- terms(object)
    Terms <- delete.response(tt)
    if (missing(newdata) || is.null(newdata)) {
        m <- model.frame(object)
        newd <- model.matrix(Terms, m)
    }
    else {
        m <- model.frame(Terms, newdata, na.action = na.action)
        newd <- model.matrix(Terms, m)
    }
    type <- match.arg(type)
    lin.pred <- as.vector(newd %*% object$coefficients)
    link <- object$link
    res <- switch(link, logit = plogis(lin.pred), probit = pnorm(lin.pred), 
        cloglog = 1 - exp(-exp(lin.pred)))
    if (type == "response") 
        pred <- res
    else pred <- lin.pred
    if (se.fit) {
        cov <- solve(object$hessian)
        var.lin.pred <- diag(newd %*% cov %*% t(newd))
        var.res <- switch(link, logit = exp(2 * lin.pred)/(1 + 
            exp(lin.pred))^4, probit = dnorm(lin.pred)^2, cloglog = (exp(-exp(lin.pred)) * 
            exp(lin.pred))^2) * var.lin.pred
        if (type == "response") 
            se <- sqrt(var.res)
        else se <- sqrt(var.lin.pred)
        if (!is.null(conf.level)) {
            alpha <- 1 - conf.level
            lower <- lin.pred - qnorm(1 - alpha/2) * sqrt(var.lin.pred)
            upper <- lin.pred + qnorm(1 - alpha/2) * sqrt(var.lin.pred)
            res.lower <- switch(link, logit = plogis(lower), 
                probit = pnorm(lower), cloglog = 1 - exp(-exp(lower)))
            res.upper <- switch(link, logit = plogis(upper), 
                probit = pnorm(upper), cloglog = 1 - exp(-exp(upper)))
            if (type == "response") {
                lwr <- res.lower
                upr <- res.upper
            }
            else {
                lwr <- lower
                upr <- upper
            }
        }
    }
    names(pred) <- 1:length(lin.pred)
    if (!is.null(conf.level)) {
        list(fit = pred, se.fit = se, lower = lwr, upper = upr)
    }
    else if (se.fit) 
        list(fit = pred, se.fit = se)
    else pred
}


residuals.gt <- function (object, type = c("deviance", "pearson", "response"), 
    ...) 
{
    type <- match.arg(type)
    r <- object$residuals
    zhat <- object$fitted.values
    z <- object$z
    res <- switch(type, response = r, pearson = r/sqrt(zhat * 
        (1 - zhat)), deviance = sqrt(-2 * (z * log(zhat) + (1 - 
        z) * log(1 - zhat))) * sign(z - zhat))
    res
}


summary.gt.mp <- function (object, ...) 
{
    coef.p <- object$coefficients
    cov.mat <- solve(object$hessian)
    dimnames(cov.mat) <- list(names(coef.p), names(coef.p))
    var.cf <- diag(cov.mat)
    s.err <- sqrt(var.cf)
    zvalue <- coef.p/s.err
    dn <- c("Estimate", "Std. Error")
    pvalue <- 2 * pnorm(-abs(zvalue))
    coef.table <- cbind(coef.p, s.err, zvalue, pvalue)
    dimnames(coef.table) <- list(names(coef.p), c(dn, "z value", 
        "Pr(>|z|)"))
    keep <- match(c("call", "link", "Gibbs.sample.size", "counts"), names(object), 0)
    ans <- c(object[keep], list(coefficients = coef.table, cov.mat = cov.mat))
    class(ans) <- "summary.gt.mp"
    return(ans)
}


print.summary.gt.mp <- function (x, digits = max(3, getOption("digits") - 3), 
       signif.stars = getOption("show.signif.stars"), 
       ...) 
{
    obj <- x
    cat("\nCall:\n")
    cat(paste(deparse(obj$call), sep = "\n", collapse = "\n"), 
        "\n", sep = "")
    cat("\nCoefficients:\n")
    coefs <- obj$coefficients
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
        na.print = "NA", ...)
    cat("\nNumber of Gibbs samples generated in each E step: ", 
        obj$Gibbs.sample.size, "\n", "Number of iterations in EM algorithm: ", 
        obj$counts, "\n", sep = "")
    cat("\n")
    invisible(obj)
}


print.gt <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    if (length(coef(x))) {
        cat("Coefficients:\n")
        print.default(format(coef(x), digits = digits), print.gap = 2, 
            quote = FALSE)
    }
    else cat("No coefficients\n")
    if (!is.null(x$df.null)) {
        cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ", 
            x$df.residual, "Residual\n")
        cat("Null Deviance:\t   ", format(signif(x$null.deviance, 
            digits)), "\nResidual Deviance:", format(signif(x$deviance, 
            digits)), "\tAIC:", format(signif(x$aic, digits)), "\n")
    }
    invisible(x)
}


gtreg.halving <- function(formula, data, groupn, subg, retest, sens = 1, spec = 1,
        linkf = c("logit", "probit", "cloglog"),
        sens.ind = NULL, spec.ind = NULL,
        start = NULL, control = gt.control(...), ...) {

    call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "groupn"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    gr <- model.extract(mf, "groupn")
    if (!is.na(pos <- match(deparse(substitute(retest)), names(data))))
        retest <- data[, pos]
    if (!is.na(pos <- match(deparse(substitute(subg)), names(data))))
        subg <- data[, pos]
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm))
            names(Y) <- nm
    }
    X <- if (!is.empty.model(mt))
        model.matrix(mt, mf)
    else matrix(, NROW(Y), 0)
    linkf <- match.arg(linkf)
    fit <-  EM.halving(Y, X, gr, subg, retest, sens, spec, linkf,
             sens.ind, spec.ind, start, control)
    fit <- c(fit, list(call = call, formula = formula, method = "Xie",
          link = linkf, terms = mt))
    class(fit) <- "gt"
    fit

}


EM.halving <- function (Y, X, groupn, subg, ret, sens, spec, linkf,
    sens.ind, spec.ind,
    start = NULL, control = gt.control())
{
    if (control$time)
        start.time <- proc.time()
    if (is.null(sens.ind))
        sens.ind <- sens
    if (is.null(spec.ind))
        spec.ind <- spec
    z <- tapply(Y, groupn, tail, n = 1)
    num.g <- max(groupn)
    K <- ncol(X)
    if (is.null(start)) {
        if (K == 1)
            cova.mean <- as.matrix(tapply(X, groupn, mean))
        else {
            temp <- by(X, groupn, colMeans)
            cova.mean <- do.call(rbind, temp)
        }
        beta.old <- lm.fit(cova.mean, z)$coefficients
    } else beta.old <- start
    sam <- length(Y)
    vec <- 1:sam
    group.sizes <- tapply(Y, groupn, length)
    diff <- 1
    counts <- 1
    extra.loop <- FALSE
    next.loop <- TRUE
    a0 <- ifelse(ret == 1, sens.ind, 1 - sens.ind)
    a1 <- ifelse(ret == 0, spec.ind, 1 - spec.ind)
    while (next.loop) {
        xib <- X %*% beta.old
        pijk <- switch(linkf, logit = plogis(xib),
            probit = pnorm(xib), cloglog = 1 -
                exp(-exp(xib)))
        erf <- 2 * pijk - 1
        prodp <- tapply(1 - pijk, groupn, prod)
        den2 <- rep(spec * prodp + (1 - sens) * (1 - prodp),
            group.sizes)
        expect <- rep(NA, times = sam)
        i <- 1
        while (i <= sam) {
            if (Y[i] == 0)
                expect[i] <- (1 - sens) * pijk[i]/den2[i]
            else {
            if (subg[i] == 0) {
                vec1 <- vec[groupn == groupn[i]]
                gs <- length(vec1)
                sub1 <- vec1[1:ceiling(gs/2)]
                sub2 <- vec1[(ceiling(gs/2) + 1):gs]
                if (subg[vec1[gs]] == 0) {
                    den <- (1-spec)*spec^2*prod(1-pijk[sub1])*prod(1-pijk[sub2])+
                           spec*(1-sens)*sens*prod(1-pijk[sub1])*(1-prod(1-pijk[sub2]))+
                           spec*(1-sens)*sens*prod(1-pijk[sub2])*(1-prod(1-pijk[sub1]))+
                           (1-sens)^2*sens*(1-prod(1-pijk[sub1]))*(1-prod(1-pijk[sub2]))
                    ab1 <- (1-sens)*sens*(spec*prod(1-pijk[sub2])+
                                  (1-sens)*(1-prod(1-pijk[sub2])))
                    ab2 <- (1-sens)*sens*(spec*prod(1-pijk[sub1])+
                                  (1-sens)*(1-prod(1-pijk[sub1])))
                    for (l1 in sub1) {
                        expect[l1]<-ab1*pijk[l1]/den
                    }
                    for (l1 in sub2) {
                        expect[l1]<-ab2*pijk[l1]/den
                    }
                }
                if (subg[vec1[gs]] == 1) {
                    mb2 <- 1
                    for (l in sub2) {
                         temp <- a0[l] * pijk[l] + a1[l] * (1 - pijk[l])
                         mb2 <- mb2 * temp
                    }
                    null <- 1
                    for (l in sub2) {
                         temp <- a1[l] * (1 - pijk[l])
                         null <- null * temp
                    }
                    den <- (1-spec)^2*spec*null*prod(1-pijk[sub1])+
                           (1-spec)*(1-sens)*sens*null*(1-prod(1-pijk[sub1]))+
                           spec*sens^2*(mb2-null)*prod(1-pijk[sub1])+
                           (1-sens)*sens^2*(mb2-null)*(1-prod(1-pijk[sub1]))
                    ab1 <- (1-sens)*sens*(mb2*sens+null*(1-sens-spec))
                    for (l1 in sub1) {
                        expect[l1]<-ab1*pijk[l1]/den
                    }
                    for (l1 in sub2) {
                        temp <- a0[l1] * pijk[l1] + a1[l1] * (1 - pijk[l1])
                        num <- mb2/temp * a0[l1] * pijk[l1] * sens^2*(spec*prod(1-pijk[sub1])+
                                  (1-sens)*(1-prod(1-pijk[sub1])))
                        expect[l1]<-num/den
                    }
                }
                i <- l1
            } else {
                vec1 <- vec[groupn == groupn[i]]
                gs <- length(vec1)
                sub1 <- vec1[1:ceiling(gs/2)]
                sub2 <- vec1[(ceiling(gs/2) + 1):gs]
                if (subg[vec1[gs]] == 0) {
                    mb2 <- 1
                    for (l in sub1) {
                         temp <- a0[l] * pijk[l] + a1[l] * (1 - pijk[l])
                         mb2 <- mb2 * temp
                    }
                    null <- 1
                    for (l in sub1) {
                         temp <- a1[l] * (1 - pijk[l])
                         null <- null * temp
                    }
                    den <- (1-spec)^2*spec*null*prod(1-pijk[sub2])+
                           (1-spec)*(1-sens)*sens*null*(1-prod(1-pijk[sub2]))+
                           spec*sens^2*(mb2-null)*prod(1-pijk[sub2])+
                           (1-sens)*sens^2*(mb2-null)*(1-prod(1-pijk[sub2]))
                    ab1 <- (1-sens)*sens*(mb2*sens+null*(1-sens-spec))
                    for (l1 in sub1) {
                        temp <- a0[l1]*pijk[l1]+a1[l1]*(1-pijk[l1])
                        num <- mb2/temp*a0[l1]*pijk[l1]*sens^2*(spec*prod(1-pijk[sub2])+
                                  (1-sens)*(1-prod(1-pijk[sub2])))
                        expect[l1]<-num/den

                    }
                    for (l1 in sub2) {
                        expect[l1]<-ab1*pijk[l1]/den
                    }
                }
                if (subg[vec1[gs]] == 1) {
                    mb2 <- 1
                    for (l in sub1) {
                         temp <- a0[l] * pijk[l] + a1[l] * (1 - pijk[l])
                         mb2 <- mb2 * temp
                    }
                    null <- 1
                    for (l in sub1) {
                         temp <- a1[l] * (1 - pijk[l])
                         null <- null * temp
                    }
                    mb2a <- 1
                    for (l in sub2) {
                         temp <- a0[l] * pijk[l] + a1[l] * (1 - pijk[l])
                         mb2a <- mb2a * temp
                    }
                    nulla <- 1
                    for (l in sub2) {
                         temp <- a1[l] * (1 - pijk[l])
                         nulla <- nulla * temp
                    }
                    den <- (1-spec)^3*null*nulla+
                           (1-spec)*sens^2*null*(mb2a-nulla)+
                           (1-spec)*sens^2*(mb2-null)*nulla+
                           sens^3*(mb2-null)*(mb2a-nulla)
                    for (l1 in sub1) {
                        temp <- a0[l1]*pijk[l1]+a1[l1]*(1-pijk[l1])
                        num <- mb2/temp*a0[l1]*pijk[l1]*sens^2*(mb2a*sens+nulla*(1-sens-spec))
                        expect[l1]<-num/den
                    }
                    for (l1 in sub2) {
                        temp <- a0[l1]*pijk[l1]+a1[l1]*(1-pijk[l1])
                        num <- mb2a/temp*a0[l1]*pijk[l1]*sens^2*(mb2*sens+null*(1-sens-spec))
                        expect[l1]<-num/den
                    }
                }
                i <- l1
            }
            }
            i <- i + 1
        }
        expect[expect > 1] <- 1
        expect[expect < 0] <- 0
        if (!extra.loop) {
            suppress <- function(w)
                if (any(grepl("non-integer #successes in a binomial glm", w)))
                   invokeRestart("muffleWarning")
            mod.fit <- withCallingHandlers(glm.fit(X, expect,
                family = binomial(link = linkf)), warning = suppress)
            diff <- max(abs((beta.old - mod.fit$coefficients)/beta.old))
            beta.old <- mod.fit$coefficients
            if (control$trace)
                cat("beta is", beta.old, "\tdiff is", diff, "\n")
            counts <- counts + 1
            if (diff <= control$tol || counts > control$maxit)
                extra.loop <- TRUE
        }
        else next.loop <- FALSE
    }
    pt1 <- switch(linkf, logit = -exp(xib)/(1 + exp(xib))^2,
        probit = sqrt(2) * xib * exp(-xib^2/2)/(sqrt(pi) * (1 -
            erf)) - 2 * exp(-xib^2)/(pi * (1 - erf)^2), cloglog = -exp(xib))
    pt2 <- switch(linkf, logit = 0, probit = (8 * exp(-xib^2/2) *
        erf + 2 * xib * sqrt(2 * pi) * erf^2 - 2 * xib * sqrt(2 *
        pi)) * exp(-xib^2/2)/((1 + erf)^2 * pi * (1 - erf)^2),
        cloglog = -(exp(xib - exp(xib)) + exp(2 * xib - exp(xib)) -
            exp(xib))/(exp(-exp(xib)) - 1)^2)
    nm <- pt1 + expect * pt2
    sign1 <- as.vector(sign(nm))
    nn <- as.vector(sqrt(abs(nm)))
    x2 <- X * nn
    m <- (t(x2) %*% (sign1 * x2))
    m1 <- 0
    i <- 1
    while (i <= sam) {
        vec1 <- vec[groupn == groupn[i]]
        gs <- length(vec1)
        if (Y[i] == 0) {
            for (j in vec1) {
                 wii <- ifelse(i == j, expect[i] - expect[i]^2, expect[i] *
                    (pijk[j] - expect[j]))
                 tim <- wii * X[i, ] %*% t(X[j, ])
                 m1 <- m1 + tim
            }
        } else {
            sub1 <- vec1[1:ceiling(gs/2)]
            sub2 <- vec1[(ceiling(gs/2) + 1):gs]
            for (i in sub1) {
                 for (j in sub1) { 
                     if (subg[j] == 0) { 
                         eii <- expect[i] * pijk[j]
                     } else { 
                         temp <- a0[j] * pijk[j] + a1[j] * (1 - pijk[j])
                         eii <- expect[i]/temp * a0[j] * pijk[j]
                     }
                     wii <- ifelse(i == j, expect[i] - expect[i]^2, eii - expect[i] * expect[j])
                     tim <- wii * X[i, ] %*% t(X[j, ])
                     m1 <- m1 + tim
                 }
                 for (j in sub2) { 
                     if (subg[j] == 0) { 
                         temp<-spec*prod(1-pijk[sub2])+(1-sens)*(1-prod(1-pijk[sub2]))
                         eii <- expect[i]*(1-sens)*pijk[j]/temp
                     } else { 
                         mb2a <- 1
                         for (l in sub2) {
                              temp <- a0[l] * pijk[l] + a1[l] * (1 - pijk[l])
                              mb2a <- mb2a * temp
                         }
                         nulla <- 1
                         for (l in sub2) {
                              temp <- a1[l] * (1 - pijk[l])
                              nulla <- nulla * temp
                         }
                         temp <- a0[j]*pijk[j]+a1[j]*(1-pijk[j])
                         tempa <- mb2a * sens + nulla * (1 - sens - spec)
                         eii <- expect[i]/tempa*sens*a0[j]*pijk[j]*mb2a/temp
                     }
                     wii <- ifelse(i == j, expect[i] - expect[i]^2, eii - expect[i] * expect[j])
                     tim <- wii * X[i, ] %*% t(X[j, ])
                     m1 <- m1 + tim
                 }
            }
            for (i in sub2) {
                 for (j in sub1) { 
                     if (subg[j] == 0) { 
                         temp<-spec*prod(1-pijk[sub1])+(1-sens)*(1-prod(1-pijk[sub1]))
                         eii <- expect[i] * (1-sens)* pijk[j]/temp
                     } else { 
                         mb2 <- 1
                         for (l in sub1) {
                              temp <- a0[l] * pijk[l] + a1[l] * (1 - pijk[l])
                              mb2 <- mb2 * temp
                         }
                         null <- 1
                         for (l in sub1) {
                              temp <- a1[l] * (1 - pijk[l])
                              null <- null * temp
                         }
                         temp <- a0[j]*pijk[j]+a1[j]*(1-pijk[j])
                         tempa <- mb2*sens+null*(1-sens-spec)
                         eii <- expect[i]/tempa*sens*a0[j]*pijk[j]*mb2/temp
                     }
                     wii <- ifelse(i == j, expect[i] - expect[i]^2, eii - expect[i] * expect[j])
                     tim <- wii * X[i, ] %*% t(X[j, ])
                     m1 <- m1 + tim
                 }
                 for (j in sub2) { 
                     if (subg[j] == 0) { 
                         eii <- expect[i] * pijk[j]
                     } else { 
                         temp <- a0[j] * pijk[j] + a1[j] * (1 - pijk[j])
                         eii <- expect[i]/temp * a0[j] * pijk[j]
                     }
                     wii <- ifelse(i == j, expect[i] - expect[i]^2, eii - expect[i] * expect[j])
                     tim <- wii * X[i, ] %*% t(X[j, ])
                     m1 <- m1 + tim
                 }
            }
        }
        i <- i + 1
    }
    H <- -(m + m1)
    zhat <- sens + (1 - sens - spec) * prodp
    residual <- z - zhat
    logl <- 0
    for (grn in 1:num.g) {
         if (z[grn] == 1) {
            vec1 <- vec[groupn == grn]
            gs <- length(vec1)
            sub1 <- vec1[1:ceiling(gs/2)]
            sub2 <- vec1[(ceiling(gs/2) + 1):gs]
            if (subg[vec1[1]] == 0) {
                if (subg[vec1[gs]] == 0) {
                    prob1 <- (1-spec)*spec^2*prod(1-pijk[sub1])*prod(1-pijk[sub2])+
                           spec*(1-sens)*sens*prod(1-pijk[sub1])*(1-prod(1-pijk[sub2]))+
                           spec*(1-sens)*sens*prod(1-pijk[sub2])*(1-prod(1-pijk[sub1]))+
                           (1-sens)^2*sens*(1-prod(1-pijk[sub1]))*(1-prod(1-pijk[sub2]))
                }
                if (subg[vec1[gs]] == 1) {
                    mb2 <- 1
                    for (l in sub2) {
                         temp <- a0[l] * pijk[l] + a1[l] * (1 - pijk[l])
                         mb2 <- mb2 * temp
                    }
                    null <- 1
                    for (l in sub2) {
                         temp <- a1[l] * (1 - pijk[l])
                         null <- null * temp
                    }
                    prob1 <- (1-spec)^2*spec*null*prod(1-pijk[sub1])+
                           (1-spec)*(1-sens)*sens*null*(1-prod(1-pijk[sub1]))+
                           spec*sens^2*(mb2-null)*prod(1-pijk[sub1])+
                           (1-sens)*sens^2*(mb2-null)*(1-prod(1-pijk[sub1]))
                }
            } else {
                if (subg[vec1[gs]] == 0) {
                    mb2 <- 1
                    for (l in sub1) {
                         temp <- a0[l] * pijk[l] + a1[l] * (1 - pijk[l])
                         mb2 <- mb2 * temp
                    }
                    null <- 1
                    for (l in sub1) {
                         temp <- a1[l] * (1 - pijk[l])
                         null <- null * temp
                    }
                    prob1 <- (1-spec)^2*spec*null*prod(1-pijk[sub2])+
                           (1-spec)*(1-sens)*sens*null*(1-prod(1-pijk[sub2]))+
                           spec*sens^2*(mb2-null)*prod(1-pijk[sub2])+
                           (1-sens)*sens^2*(mb2-null)*(1-prod(1-pijk[sub2]))
                }
                if (subg[vec1[gs]] == 1) {
                    mb2 <- 1
                    for (l in sub1) {
                         temp <- a0[l] * pijk[l] + a1[l] * (1 - pijk[l])
                         mb2 <- mb2 * temp
                    }
                    null <- 1
                    for (l in sub1) {
                         temp <- a1[l] * (1 - pijk[l])
                         null <- null * temp
                    }
                    mb2a <- 1
                    for (l in sub2) {
                         temp <- a0[l] * pijk[l] + a1[l] * (1 - pijk[l])
                         mb2a <- mb2a * temp
                    }
                    nulla <- 1
                    for (l in sub2) {
                         temp <- a1[l] * (1 - pijk[l])
                         nulla <- nulla * temp
                    }
                    prob1 <- (1-spec)^3*null*nulla+
                           (1-spec)*sens^2*null*(mb2a-nulla)+
                           (1-spec)*sens^2*(mb2-null)*nulla+
                           sens^3*(mb2-null)*(mb2a-nulla)
                }
            }
         } else prob1 <- 1 - zhat[grn]
         logl <- logl - log(prob1)
    }
    aic <- 2 * logl + 2 * K
    if (diff > control$tol && counts > control$maxit)
        warning("EM algorithm did not converge.")
    if (control$time) {
        end.time <- proc.time()
        save.time <- end.time - start.time
        cat("\n Number of minutes running:", round(save.time[3]/60, 2), "\n \n")
    }
    list(coefficients = beta.old, hessian = H, fitted.values = zhat,
        deviance = 2 * logl, aic = aic,
        counts = counts - 1, residuals = residual, z = z)
}


sim.halving <- function (x = NULL, gshape = 20, gscale = 2, par,
    linkf = c("logit", "probit", "cloglog"),
    sample.size, group.size, sens = 1, spec = 1,
    sens.ind = NULL, spec.ind = NULL)
{
    if (is.null(sens.ind))
        sens.ind <- sens
    if (is.null(spec.ind))
        spec.ind <- spec
    if (is.null(x)) {
        x <- rgamma(n = sample.size, shape = gshape, scale = gscale)
        X <- cbind(1, x)
    }
    else {
        X <- cbind(1, x)
        sample.size <- nrow(X)
    }
    linkf <- match.arg(linkf)
    pijk <- switch(linkf, logit = plogis(X %*% par),
            probit = pnorm(X %*% par),
            cloglog = 1 - exp(-exp(X %*% par)))
    ind <- rbinom(n = sample.size, size = 1, prob = pijk)
    num.g <- ceiling(sample.size/group.size)
    vec <- 1:sample.size
    groupn <- rep(1:num.g, each = group.size)[vec]
    save.sum <- tapply(ind, groupn, sum)
    save.group <- as.vector(ifelse(save.sum > 0, 1, 0))
    save.obs <- rep(NA, num.g)
    subgroup <- ret <- rep(NA, sample.size)
    for (grn in 1:num.g) {
         vec1 <- vec[groupn == grn]
         gs <- length(vec1)
         save.obs[grn] <- ifelse(save.group[grn] == 1, rbinom(1, 1, sens),
                1 - rbinom(1, 1, spec))
         if (save.obs[grn] == 1) {
             sub1 <- vec1[1:ceiling(gs/2)]
             sub2 <- vec1[(ceiling(gs/2) + 1):gs]
             tZ1 <- sum(ind[sub1])
             tZ2 <- sum(ind[sub2])
             Z1 <- ifelse(tZ1 == 1, rbinom(1,
                 1, sens), 1 - rbinom(1, 1, spec))
             Z2 <- ifelse(tZ2 == 1, rbinom(1,
                 1, sens), 1 - rbinom(1, 1, spec))
             if (Z1 == 1) {
                 for (i1 in sub1) {
                     ret[i1] <- ifelse(ind[i1] == 1, rbinom(1,
                        1, sens.ind), 1 - rbinom(1, 1, spec.ind))
                 }
             }
             if (Z2 == 1) {
                 for (i1 in sub2) {
                     ret[i1] <- ifelse(ind[i1] == 1, rbinom(1,
                        1, sens.ind), 1 - rbinom(1, 1, spec.ind))
                 }
             }
             subgroup[sub1] <- Z1
             subgroup[sub2] <- Z2
         }
    }
    gres <- rep(save.obs, each = group.size)[vec]
    grd <- data.frame(gres = gres, x = x, groupn = groupn, ind = ind,
            retest = ret, subgroup = subgroup)
    if (ncol(X) > 2)
        for (i in 2:ncol(X))
             colnames(grd)[i] <- paste("x", i - 1, sep="")
    grd
}
