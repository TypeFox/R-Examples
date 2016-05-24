bruto <-
function (x, y, w = rep(1, n), wp = rep(1/np, np), dfmax, cost = 2, 
          maxit.select = 20, maxit.backfit = 20, thresh = 1e-04,
          trace.bruto = FALSE, start.linear = TRUE, fit.object, ...)
{
    this.call <- match.call()
    y <- as.matrix(y)
    x <- as.matrix(x)
    np <- ncol(y)
    d <- dim(x)
    nq <- d[2]
    n <- d[1]
    xnames <- dimnames(x)[[2]]
    if (!length(xnames)) 
        xnames <- NULL
    ynames <- dimnames(y)[[2]]
    if (!length(ynames)) 
        ynames <- NULL
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"
    storage.mode(w) <- "double"
    storage.mode(wp) <- "double"
    storage.mode(cost) <- "double"
    if (missing(fit.object)) {
        nknotl <- function(n) {
            a1 <- log(50)/log(2)
            a2 <- log(100)/log(2)
            a3 <- log(140)/log(2)
            a4 <- log(200)/log(2)
            cx <- as.numeric(cut(n, c(0, 50, 200, 800, 3200)))
            if (is.na(cx)) 
                cx <- 5
            floor(switch(cx, n, 2^(a1 + ((a2 - a1) * (n - 50))/150), 
                2^(a2 + ((a3 - a2) * (n - 200))/600), 2^(a3 + 
                  ((a4 - a3) * (n - 800))/2400), 200 + (n - 3200)^0.2) + 
                6)
        }
        check.range <- apply(x, 2, var)
        if (any(check.range < .Machine$double.eps)) 
            stop(paste("A column of x is constant;",
                       "do not include an intercept column"))
        nkmax <- nknotl(n) - 4
        coef <- matrix(double(nkmax * np * nq), ncol = nq)
        ybar <- apply(y * w, 2, sum)/sum(w)
        if (start.linear && (nq * cost > n)) 
            start.linear <- FALSE
        if (start.linear) {
            start.fit <- polyreg(x, y, w)
            eta <- fitted(start.fit)
            coef[seq(from = 2, by = 2, length = np), ] <-
                t(start.fit$coef)[, -1]
            type <- as.integer(rep(2, nq))
            df <- as.double(rep(1, nq))
        }
        else {
            eta <- outer(rep(1, n), ybar)
            type <- integer(nq)
            df <- double(nq)
        }
        nk <- integer(nq)
        knot <- matrix(double((nkmax + 4) * nq), ncol = nq)
        Match <- matrix(integer(n * nq), ncol = nq)
        nef <- integer(nq)
        lambda <- double(nq)
        xrange <- matrix(double(2 * nq), 2, nq)
        df <- double(nq)
        if (missing(dfmax)) 
            dfmax <- (2 * nkmax)/3
        if (length(dfmax) != nq) 
            dfmax <- rep(dfmax, length = nq)
        if (cost > 0) {
            TD <- (n - sum(df))/cost
            TT <- dfmax > TD
            if (any(TT)) 
                dfmax[TT] <- TD
        }
        storage.mode(dfmax) <- "double"
    }
    else {
        this.call <- fit.object$call
        ybar <- fit.object$ybar
        nkmax <- fit.object$nkmax
        dfmax <- fit.object$dfmax
        eta <- fit.object$fitted.values
        if (is.null(eta)) 
            eta <- predict(fit.object, x)
        nk <- fit.object$nk
        knot <- fit.object$knot
        Match <- fit.object$Match
        nef <- fit.object$nef
        lambda <- fit.object$lambda
        coef <- fit.object$coef
        type <- unclass(fit.object$type)
        xrange <- fit.object$xrange
        maxit.select <- 0
        maxit.backfit <- fit.object$nit[2]
        df <- fit.object$df
    }
    maxit <- as.integer(c(maxit.select, maxit.backfit))
    names(df) <- xnames
    names(maxit) <- c("selection", "backfitting")
    gcv.select <- if (maxit.select) 
        matrix(double(maxit.select * nq), nq, maxit.select)
    else double(1)
    gcv.backfit <- if (maxit.backfit) 
        matrix(double(maxit.backfit * nq), nq, maxit.backfit)
    else double(1)
    df.select <- if (maxit.select) 
        matrix(double(maxit.select * nq), nq, maxit.select)
    else double(1)
    names(lambda) <- xnames
    fit <-
        .Fortran("bruto",
                 x,
                 as.integer(n),
                 as.integer(nq), 
                 y,
                 as.integer(np),
                 w,
                 knot = knot,
                 nkmax = as.integer(nkmax), 
                 nk = nk,
                 wp,
                 Match = Match,
                 nef = nef,
                 dfmax = dfmax, 
                 cost = cost,
                 lambda = lambda,
                 df = df,
                 coef = coef,
                 type = type, 
                 xrange = xrange,
                 gcv.select = gcv.select,
                 gcv.backfit = gcv.backfit, 
                 df.select = df.select,
                 maxit = maxit,
                 nit = maxit,
                 fitted.values = eta, 
                 residuals = y - eta,
                 as.double(thresh),
                 double((2 * np + 2) * ((n + 1) + 1) + (2 * np + 16) *
                        (n + 1) + 2 * (n + 1) + np),
                 integer(n),
                 trace.bruto,
                 PACKAGE = "mda")[c("knot", "nkmax", "nk", "Match",
                 "nef", "dfmax", "cost", "lambda", "df", "coef", "type",
                 "xrange", "gcv.select", "gcv.backfit", "df.select",
                 "maxit", "nit", "fitted.values", "residuals")]
    if (TN <- fit$nit[1]) {
        TT <- fit$gcv.select[, seq(TN), drop = FALSE]
        dimnames(TT) <- list(xnames, NULL)
    }
    else TT <- NULL
    fit$gcv.select <- TT
    if (TN) {
        TT <- fit$df.select[, seq(TN), drop = FALSE]
        dimnames(TT) <- list(xnames, NULL)
    }
    else TT <- NULL
    fit$df.select <- TT
    if (TN <- fit$nit[2]) {
        TT <- fit$gcv.backfit[, seq(TN), drop = FALSE]
        dimnames(TT) <- list(xnames, NULL)
    }
    else TT <- NULL
    fit$gcv.backfit <- TT
    TT <- factor(fit$type, levels = 1:3, labels = c("excluded", 
        "linear", "smooth"))
    names(TT) <- xnames
    fit$type <- TT
    fit$ybar <- ybar
    fit$call <- this.call
    structure(fit, class = "bruto")
}

