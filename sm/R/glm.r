"sm.binomial" <- function (x, y, N = rep(1, length(y)), h, ...) {
    x.name <- deparse(substitute(x))
    y.name <- deparse(substitute(y))
    opt <- sm.options(list(...))
    if (any(is.na(c(x, y, N)))) {
        xy <- cbind(x, y, N)
        ok <- as.logical(apply(!is.na(xy), 1, prod))
        xy <- xy[ok, ]
        x <- as.vector(xy[, 1])
        y <- as.vector(xy[, 2])
        N <- as.vector(xy[, 3])
        if(opt$verbose>0) cat("warning: missing data are removed\n")
    }
    n <- length(y)
    replace.na(opt, display, "line")
    replace.na(opt, ngrid, 25)
    replace.na(opt, ylim, c(0, 1))
    replace.na(opt, col,  "black")
    replace.na(opt, nbins, round((n > 100) * 8 * log(n)))
    display <- opt$display
    if (length(x) != n)
        stop("x and y have different length")
    if (length(N) != n)
        stop("N and y have different length")
    y <- as.integer(y)
    if (min(diff(x)) < 0) {
        y <- y[order(x)]
        x <- sort(x)
    }
    replace.na(opt, eval.points, seq(min(x), max(x), length = opt$ngrid))
    replace.na(opt, nbins, round((nobs > 100) * 8 * log(n)/ndim))
    if (all(N == 1))
        yplot <- jitter(y, amount = 0)
    else yplot <- y/N
    if (display != "none" & opt$add == FALSE) {
        replace.na(opt, xlab, x.name)
        replace.na(opt, ylab, paste("Pr{", y.name, "}", sep = ""))
        plot(x, yplot, ylim = opt$ylim, xlab = opt$xlab, ylab = opt$ylab,
            col = 1, type = "n")
        abline(0, 0, col = 1, lty = 3)
        abline(1, 0, col = 1, lty = 3)
    }
    if (display != "none")
        points(x, yplot, pch = opt$pch, col = opt$col)
    rawdata <- list(x = x, y = y, N = N, nbins = opt$nbins, nobs = n, ndim = 1)
    if (opt$nbins > 0) {
        bins <- binning(x, y, nbins = opt$nbins)
        binsN <- binning(x, N, nbins = opt$nbins)
        x <- bins$x
        y <- round(bins$sums)
        N <- round(binsN$sums)
        nx <- length(y)
    }
    result <- sm.glm(x, cbind(y, N - y), family = binomial(),
        h = h, eval.points = opt$eval.points, start = log((y +
            0.5)/(N - y + 1)), options=opt)
    result$call <- match.call()
    if (display != "none") {
        lines(result$eval.points, result$estimate, col = opt$col)
        if (display == "se") {
            lines(result$eval.points, result$lower, lty = 3, col = opt$col)
            lines(result$eval.points, result$upper, lty = 3, col = opt$col)
        }
    }
    result$data <- list(x = x, y = y, N = N, nbins = opt$nbins)
    invisible(result)
}

"sm.binomial.bootstrap" <- function (x, y, N = rep(1, length(x)), h, degree = 1,
          fixed.disp = FALSE, ...) {
    rbetabinom <- function(n, size, prob, disp) {
        if (disp > 1 & min(size) > 2 & min(size) > disp) {
            psi <- (disp - 1)/(size - 1)
            alpha <- prob * (1/psi - 1)
            beta <- (1 - prob) * (1/psi - 1)
            p <- rbeta(n, alpha, beta)
            y <- rbinom(n, size, p)
        }
        else y <- rbinom(n, size, prob)
        return(y)
    }
    family<- binomial()
    opt <- sm.options(list(...))
    verbose <- opt$verbose
    nboot <- opt$nboot
    D <- function (mu, y, wt) sum(family$dev.resids(y, mu, wt))
    n <- length(x)
    sm <- sm.binomial(x, y, N, h, xlab = deparse(substitute(x)),
        ylab = paste("Pr{", deparse(substitute(y)), "}", sep = ""), ...)
    X <- cbind(1, poly(x, degree))
    colnames(X) <- seq(len=ncol(X))
    glm.model <- glm.fit(X, cbind(y, N - y), family = family)
    glm.fitted <- fitted(glm.model)
    glm.resid <- residuals(glm.model)
    lines(x, glm.fitted, lty = 2, col = 2)
    p.boot <- 0
    sm.orig <- sm.binomial(x, y, N, h, eval.points = x, display = "none")
    sm.fitted <- sm.orig$estimate
    disp.orig <- D(sm.fitted, y/N, N)/(n - degree - 1)
    if (fixed.disp) disp <- 1
    else disp <- disp.orig
    ts.orig <- (D(glm.fitted, y/N, N) - D(sm.fitted, y/N, N))/disp
    if(verbose>0) {
      cat("Dispersion parameter = ", disp.orig, "\n")
      cat("Test statistic = ", ts.orig, "\n")
    }
    yboot <- rep(NA, n)
    if(verbose>1) cat("Running to ",nboot,": ")
    for (i in 1:nboot) {
        yboot <- rbetabinom(n, N, glm.fitted, disp)
        if(verbose>1) cat("Sample:", i, " ")
        sm.fit <- sm.glm(x, cbind(yboot, N - yboot), family = family,
            h, eval.points = x, start = log((yboot + 0.5)/(N - yboot + 0.5)))
        sm.fitted <- sm.fit$estimate
        ts.boot <- (D(glm.fitted, yboot/N, N) - D(sm.fitted, yboot/N, N))/disp
        if (ts.boot > ts.orig) p.boot <- p.boot + 1
        lines(x, sm.fitted, lty = 2, col = 4)
    }
    if(verbose>1) cat("\n")
    lines(sm$eval.points, sm$estimate)
    p.boot <- p.boot/(nboot + 1)
    if(verbose>0) cat("Observed significance = ", p.boot, "\n")
    invisible(list(call = match.call(), significance = p.boot,
        test.statistic = ts.orig, dispersion = disp.orig))
}

"sm.poisson" <- function (x, y, h, ...) {
    x.name <- deparse(substitute(x))
    y.name <- deparse(substitute(y))
    opt <- sm.options(list(...))
    verbose <- opt$verbose    
    if (any(is.na(c(x, y)))) {
        xy <- cbind(x, y)
        ok <- as.logical(apply(!is.na(xy), 1, prod))
        xy <- xy[ok, ]
        y <- as.vector(xy[, ncol(xy)])
        x <- xy[, -ncol(xy), drop = TRUE]
        if(opt$verbose>0)
          cat("warning: missing data are removed\n")
    }
    y <- as.integer(y)
    n <- length(y)
    replace.na(opt, display, "line")
    replace.na(opt, ngrid, 25)
    replace.na(opt, ylim, c(0, 1))
    replace.na(opt, pch, 1)
    replace.na(opt, col, 2)
    replace.na(opt, nbins, round((n > 100) * 8 * log(n)))
    display <- opt$display
    if (min(diff(x)) < 0) {
        y <- y[order(x)]
        x <- sort(x)
    }
    if (!(opt$display %in% "none") & opt$add %in% FALSE) {
        replace.na(opt, xlab, x.name)
        replace.na(opt, ylab, y.name)
        plot(x, y, xlab = opt$xlab, ylab = opt$ylab, col = 1, type = "n")
    }
    if (display != "none")
        points(x, y, pch = opt$pch, col = opt$col)
    replace.na(opt, eval.points, seq(min(x), max(x), length = opt$ngrid))
    rawdata <- list(x = x, y = y, nbins = opt$nbins, nobs = n,
        ndim = 1)
    if (opt$nbins > 0) {
        bins <- binning(x, y, nbins = opt$nbins)
        x <- bins$x
        y <- round(bins$sums)
        nx <- length(y)
        freq <- bins$x.freq
    }
    else freq <- rep(1, n)
    result <- sm.glm(x, y, family = poisson(), h = h,
                     eval.points = opt$eval.points, start = log(pmax(0.167, y)),
                     offset = log(freq), options=opt)
    result$call <- match.call()
    if (display != "none") {
        lines(result$eval.points, result$estimate, col = opt$col)
        if (display == "se") {
            lines(result$eval.points, result$lower, lty = 3, col = opt$col)
            lines(result$eval.points, result$upper, lty = 3, col = opt$col)
        }
    }
    result$data <- list(x = x, y = y, weights = freq, nbins = opt$nbins)
    invisible(result)
}

"sm.poisson.bootstrap" <- function (x, y, h,  degree = 1, fixed.disp = FALSE, intercept = TRUE, ...) {
    rNegBin <- function(n, mean, disp) {
        if (disp > 1) {
            p <- 1/disp
            r <- mean/(disp - 1)
            theta <- (rgamma(n, r) * (1 - p))/p
            y <- rpois(n, theta)
        }
        else y <- rpois(n, mean)
        return(y)
    }
    family<- poisson()
    opt <- sm.options(list(...))
    verbose <- opt$verbose
    nboot <- opt$nboot
    D <- function(mu, y, w, residuals = FALSE)
           sum(family$dev.resids(y, mu, w))
    x.name <- deparse(substitute(x))
    y.name <- deparse(substitute(y))
    y <- as.integer(y)
    y <- y[order(x)]
    x <- sort(x)
    sm <- sm.poisson(x, y, h, xlab = x.name, ylab = y.name, col = 3, ...)
    if (intercept)
        X <- cbind(1, poly(x, degree))
    else X <- outer(x, 1:degree, "^")
    colnames(X) <- seq(len=ncol(X))
    glm.model <- glm.fit(X, y, family = family)
    glm.fitted <- fitted(glm.model)
    lines(x, glm.fitted, col = 5)
    p.boot <- 0
    sm.orig <- sm.poisson(x, y, h, eval.points = x, display = "none", ...)
    sm.fitted <- sm.orig$estimate
    disp.orig <- D(sm.fitted, y, 1)/(length(y) - ncol(X))
    if (fixed.disp)
        disp <- 1
    else disp <- disp.orig
    ts.orig <- (D(glm.fitted, y, 1) - D(sm.fitted, y, 1))/disp
    if(verbose>0){
      cat("Dipersion parameter = ", disp.orig, "\n")
      cat("Test statistic = ", ts.orig, "\n")
    }
    if(verbose>1) cat("Running to ",nboot,": ")
    for (i in 1:nboot) {
        if(verbose>1) cat(i, " ")
        yboot <- rNegBin(length(glm.fitted), glm.fitted, disp)
        sm <- sm.poisson(x, yboot, h, eval.points = x, display = "none")
        sm.fitted <- sm$estimate
        ts.boot <- (D(glm.fitted, yboot, 1) - D(sm.fitted, yboot, 1))/disp
        if (ts.boot > ts.orig)
            p.boot <- p.boot + 1
        lines(x, sm.fitted, lty = 2, col = 6)
    }
    if(verbose>1) cat("\n")
    lines(sm$eval.points, sm$estimate, col = 3)
    lines(x, glm.fitted, col = 5)
    p.boot <- p.boot/(nboot + 1)
    if(verbose>0) cat("Observed significance = ", p.boot, "\n")
    invisible(list(call = match.call(), test.statistic = ts.orig,
        significance = p.boot, disp = disp.orig))
}

"sm.glm" <- function (x, y, family, h, eval.points, start, offset, options=list()) {
    opt <- sm.options(options)
    n <- length(x)
    verbose<- as.integer(opt$verbose)    
    X <- cbind(rep(1, n + 1), c(x, 0))
    ## in R, avoid zero weight
    if (isMatrix(y)) Y <- rbind(y, rep(1, ncol(y)))
    else Y <- c(y, 0)
    start <- c(start, 0)
    neval <- length(eval.points)
    if (missing(offset)) offset <- rep(0, n)
    W <- matrix(rep(eval.points, rep(n, neval)), ncol = n, byrow = TRUE)
    W <- W - matrix(rep(X[1:n, 2], neval), ncol = n, byrow = TRUE)
    W <- exp(-0.5 * (W/h)^2)
    if(verbose>1) cat("Cycles per point: ")
    est <- LP <- st.err <- dev <- var.eta <- rep(NA, neval)
    for (k in 1:neval) {
        X[n + 1, 2] <- eval.points[k]
        colnames(X) <- 1:2
        ## weight 0 does not work in R
        fit <- glm.fit(X, Y, weights = c(W[k, ], 1e-8), family = family,
            etastart = start, offset = c(offset, 0))
        start <- fit$linear.predictors
        LP[k] <- start[n + 1]
        dev[k] <- fit$deviance
        if(verbose>1) cat(fit$iter, " ")
        s <- W[k, ]
        mu <- fit$fitted.values[1:n]
        Wk <- diag(s * fit$weights[1:n])
        XXinv <- solve(t(X[1:n, ]) %*% Wk %*% X[1:n, ])
        Li <- XXinv %*% t(X[1:n, ]) %*% diag(s)
        var.Bi <- Li %*% diag(fit$weights[1:n]) %*% t(Li)
        var.eta[k] <- t(X[n + 1, ]) %*% var.Bi %*% as.vector(X[n + 1, ])
    }
    if(verbose>1) cat("\n")
    st.err <- sqrt(var.eta)
    est <- family$linkinv(LP)
    result <- list(call = match.call(), eval.points = eval.points,
        estimate = est, lower = family$linkinv(LP - 2 * st.err),
        upper = family$linkinv(LP + 2 * st.err), linear.predictor = LP,
        se = st.err, deviance = dev)
    invisible(result)
}

"sm.autoregression" <- function (x, h = hnorm(x), d = 1, maxlag = d, lags, se = FALSE, ask = TRUE) {
    sm.autoregression.1d <- function(x, h, x.name, lags, se = FALSE, ask = FALSE) {
        n <- length(x)
        if (any(diff(lags)) < 0)
            stop("lags must be in increasing order")
        x2.name <- paste(x.name, "(t)", sep = "")
        xlow <- min(x) - diff(range(x))/20
        xhi <- max(x) + diff(range(x))/20
        lags <- sort(lags)
        for (m in lags) {
            x1 <- x[(m + 1):n]
            x0 <- x[1:(n - m)]
            r <- sm.regression.eval.1d(x0, x1, h = h, model = "none",
                options = list(hmult = 1))
            x1.name <- paste(x.name, "(t-", as.character(m),
                ")", sep = "")
            plot(x0, x1, xlim = c(xlow, xhi), ylim = c(xlow,
                xhi), xlab = x1.name, ylab = x2.name)
            lines(r$eval.points, r$estimate)
            if (se) {
                rho1 <- acf(x0, lag.max = 1, plot = FALSE)$acf[2]
                lines(r$eval.points, r$estimate + 2 * r$se/sqrt(1 - rho1),
                      lty = 3)
                lines(r$eval.points, r$estimate - 2 * r$se/sqrt(1 - rho1),
                      lty = 3)
            }
            title(paste("Regression of ", x.name, " on past data",
                        sep = ""))
            if (ask & (m < lags[length(lags)]))
                pause()
        }
        invisible(r)
    }
    sm.autoregression.2d <- function(x, h, x.name, lags, ask = ask,
        ngrid = 20, display = "none") {
        if (dim(lags)[2] != 2) stop("dim(lags)[2] must be 2")
        evpt <- seq(quantile(x, 0.1), quantile(x, 0.9), length = ngrid)
        n <- length(x)
        nplot <- dim(lags)[1]
        for (ip in 1:nplot) {
            m1 <- min(lags[ip, ])
            m2 <- max(lags[ip, ])
            x0 <- x[1:(n - m2)]
            x1 <- x[(m2 - m1 + 1):(n - m1)]
            x2 <- x[(m2 + 1):n]
            r <- sm.regression.eval.2d(cbind(x0, x1), x2, h = c(h,
                h), model = "none", eval.points = cbind(evpt,
                evpt), weights = rep(1, n - m2), options = list(hmult = 1,
                h.weights = rep(1, n - m2), poly.index = 1))
            persp(evpt, evpt, r)
            head <- paste("Regression of ", x.name, " on past data (lags: ",
                as.character(m1), ", ", as.character(m2), ")",
                sep = "")
            title(head)
            if (ask & (ip < nplot)) pause()
        }
        invisible(r)
    }
    x.name <- deparse(substitute(x))
    if (missing(lags)) {
       if (d == 1)
          lags <- (1:maxlag)
       else 
          lags <- cbind(1:(maxlag - 1), 2:maxlag)
       }
    else
       if (isMatrix(lags)) d <- 2
    x <- as.vector(x)
    if (d == 1)
       r <- sm.autoregression.1d(x, h, x.name, lags, se = se, ask = ask)
    else 
       r <- sm.autoregression.2d(x, h, x.name, lags, ask = ask)
    invisible(r)
}


"sm.regression.autocor" <- function (x = 1:n, y, h.first, 
                      minh, maxh, method = "direct", ...) {
    GCV <- function(h, x, y, R, sqrt.R) {
        W <- sm.weight(x, x, h, options = list(hmult = 1))
        r <- (y - W %*% as.matrix(y))
        rss <- sum(r^2)
        Trace <- sum(diag(W))
        gcv.0 <- rss/(1 - Trace/length(x))^2
        Trace <- sum(diag(W %*% R))
        gcv.r <- rss/(1 - Trace/length(x))^2
        rw <- backsolve(sqrt.R, r)
        Trace <- sum(diag(W))
        gcv.ri <- sum(rw^2)/(1 - Trace/length(x))^2
        c(gcv.0, gcv.r, gcv.ri)
    }
    opt <- sm.options(list(...))
    verbose <- as.integer(opt$verbose)
    replace.na(opt, display, "plot")
    replace.na(opt, ngrid, 15)
    ngrid <- opt$ngrid
    n <- length(y)
    if (length(x) != n)
        stop("x and y must have equal length\n")
    if (missing(minh) & missing(x))
        minh <- 0.5
    if (missing(maxh) & missing(x))
        maxh <- 10
    w <- sm.weight(x, x, h = h.first, options = list(hmult = 1))
    ym <- as.vector(w %*% y)
    r <- (y - ym)
    autocov <- rep(0, n)
    for (k in 0:2) {
        u <- r[1:(n - k)] * r[(k + 1):n]
        autocov[k + 1] <- sum(u)/n
    }
    var <- autocov[1]
    rho1 <- autocov[2]/var
    rho2 <- autocov[3]/var
    a1 <- rho1 * (1 - rho2)/(1 - rho1^2)
    a2 <- (rho2 - rho1^2)/(1 - rho1^2)
    if(verbose>0) cat("AR[1:2] coeff: ", c(a1, a2), "\n")
    for (k in 3:(n - 1)) autocov[k + 1] <- a1 * autocov[k] +
        a2 * autocov[k - 1]
    autocorr <- autocov/var
    R <- diag(n)
    R <- outer(1:n, 1:n, function(i, j, r) r[abs(i - j) + 1],
        r = autocorr)
    sqrt.R <- chol(R)
    hvector <- seq(minh, maxh, length = ngrid)
    min.gcv <- Inf
    h.opt <- 0
    result <- matrix(0, ngrid, 3, dimnames = list(NULL, c("no.cor",
        "direct", "indirect")))
    if(verbose>1)
      cat(paste("Search for h (runs up to ", as.character(ngrid),
          "): ", sep = "", collapse = NULL))
    for (i in 1:ngrid) {
        h <- hvector[i]
        result[i, ] <- GCV(h, x, y, R, sqrt.R)
        cat(" ")
        cat(i)
    }
    if(verbose>1) cat("\n")
    if (!(opt$display %in% "none")) {
        maxlag <- min(30, n - 1)
        acf <- array(autocorr[1:(maxlag + 1)], dim = c(maxlag +
            1, 1, 1))
        lag <- array(0:maxlag, dim = c(maxlag + 1, 1, 1))
#        acf.plot(list(acf = acf, lag = lag, type = "correlation",
#            series = "residuals from preliminary smoothing", n.used = n))
        plot(lag, acf, sub="residuals from preliminary smoothing", type="h")
        pause()
        plot(c(hvector[1], hvector[ngrid]), c(min(result), max(result)),
            type = "n", xlab = "h", ylab = "Generalised cross-validation")
        title(paste("GCV criterion, method:", method, collapse = NULL))
        lines(hvector, result[, method], col = 2)
        pause()
    }
    h1 <- hvector[order(result[, method])[1]]
    if(verbose>0) cat("Suggested value of h: ", h1, "\n")
    sm1 <- sm.regression.eval.1d(x, y, h = h1, model = "none",
        options = list(hmult = 1))
    if (missing(x))
        x.name <- "time"
    else x.name <- deparse(substitute(x))
    if (opt$display != "none") {
        plot(x, y, xlab = x.name, ylab = deparse(substitute(y)),
            ...)
        lines(sm1$eval.points, sm1$estimate, col = 2)
    }
    sm1$aux <- list(h.first = h.first, first.sm = ym, acf = autocorr,
        raw.residuals = r)
    invisible(sm1)
}


"sm.rm" <- function (Time, y, minh = 0.1, maxh = 2, optimize = FALSE,
          rice.display = FALSE, ...) {
    rice <- function(h, nSubj, Time, ym, var, r, poly.index = 1) {
        nTime <- length(Time)
        w <- sm.weight(Time, Time, h, options = list(poly.index = poly.index))
        fitted <- w %*% ym
        rss <- sum((ym - fitted)^2)
        Trace <- sum(diag(w %*% r))
        criterion <- sqrt(rss/nTime - (var/nSubj) * (1 - 2 *
            Trace/nTime))
        criterion
    }
    if (!isMatrix(y))
        stop("y must be a matrix")
    opt <- sm.options(list(...))
    verbose <- as.integer(opt$verbose)
    replace.na(opt, ngrid, 20)
    ngrid <- opt$ngrid
    nSubj <- dim(y)[1]
    nTime <- dim(y)[2]
    if (missing(Time)) Time <- 1:nTime
    ym <- apply(y, 2, mean)
    z <- y - matrix(ym, nrow = nSubj, ncol = nTime, byrow = TRUE)
    autocov <- rep(0, nTime)
    for (k in 0:(nTime - 1)) {
        u <- z[, 1:(nTime - k)] * z[, (k + 1):nTime]
        autocov[k + 1] <- sum(u)/(nSubj * nTime)
    }
    var <- autocov[1]
    autocorr <- autocov/var
    if(verbose>0) {
       cat("Autocovariances & autocorrelations:\n")
       print(matrix(cbind(autocov, autocorr), ncol = 2,
             dimnames = list(0:(nTime - 1), c("auto-cov", "auto-corr"))))
      }                    
    r <- diag(nTime)
    for (k in 1:nTime) {
        for (j in 1:nTime) r[k, j] <- autocorr[abs(k - j) + 1]
    }
    hvector <- seq(minh, maxh, length = ngrid)
    min.obj <- Inf
    h.opt <- 0
    if(verbose>0) {
      cat("       Rice's criterion:\n")
      cat("       h    indept.   depend.\n")
    }
    result <- matrix(0, ngrid, 2, dimnames = list(NULL, c("indept",
        "depend")))
    for (i in 1:ngrid) {
        h <- hvector[i]
        obj.0 <- rice(h, nSubj, Time, ym, var, diag(nTime), opt$poly.index)
        obj.r <- rice(h, nSubj, Time, ym, var, r, opt$poly.index)
        result[i, 1] <- obj.0
        result[i, 2] <- obj.r
        if (obj.r < min.obj) {
            min.obj <- obj.r
            h.opt <- h
        }
        if(verbose>0) print(c(h, obj.0, obj.r))
    }
    if (rice.display) {
        plot(c(hvector[1], hvector[ngrid]), c(min(result), max(result)),
            type = "n", xlab = "h", ylab = "sqrt(rice criterion)")
        title(main = "Modified Rice criterion for selecting h",
            sub = paste("dashed line: assume independence,",
                " continuous: allow for correlation", collapse = NULL))
        lines(hvector, result[, 1], lty = 3)
        lines(hvector, result[, 2], lty = 1)
        pause()
    }
    if (optimize) {
        if(verbose>0)  cat("Search for optimum h using optim...\n")
        optimum <- optim(par = h.opt, fn = rice, method = "L-BFGS",
            lower = 0, nSubj = nSubj, Time = Time, ym = ym, var = var, r = r)
        print(optimum$par)
        h.opt <- optimum$par
    }
    if(verbose>0) cat("h: ", h.opt, "\n")
    if (opt$display %in% "se")
        display1 <- "line"
    else display1 <- opt$display
    sm <- sm.regression(Time, ym, h = h.opt, hmult = 1, display = display1,
            ylab = paste(deparse(substitute(y)), "(mean values)",
            collapse = NULL), add = opt$add)
    if (opt$display %in% "se") {
        W <- sm.weight(Time, sm$eval.points, h = h.opt, options = list())
        V <- (var/nSubj) * r
        se <- sqrt(diag(W %*% V %*% t(W)))
        lines(sm$eval.points, sm$estimate + 2 * se, lty = 3)
        lines(sm$eval.points, sm$estimate - 2 * se, lty = 3)
    }
    sm$aux <- list(mean = ym, var = var, autocorr = autocorr, h = h.opt)
    invisible(sm)
}


"sm.ts.pdf" <- function (x, h = hnorm(x), lags, maxlag = 1, ask = TRUE) {
    if (missing(lags))
        lags <- (1:maxlag)
    else maxlag <- max(lags)
    if (any(diff(lags) < 0))
        stop("lags must be in increasing order")
    x.name <- deparse(substitute(x))
    x <- as.vector(x)
    n <- length(x)
    marginal <- sm.density(x, ylab = "Marginal density", xlab = x.name)
    if (ask)
        pause()
    for (m in lags) {
        x1 <- x[(m + 1):n]
        x0 <- x[1:(n - m)]
        biv <- sm.density(cbind(x0, x1), h = rep(h, 2), xlab = paste(x.name,
            "(t-", as.character(m), ")", sep = ""), ylab = paste(x.name,
            "(t)", sep = ""))
        biv$lag <- m
        title(paste("Density of lagged data of ", x.name, " (lag=",
            as.character(m), ")", sep = ""))
        if (ask & (m < maxlag))
            pause()
    }
    invisible(list(marginal = marginal, bivariate = biv))
}

