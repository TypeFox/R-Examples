biexp <- function (conc, time, log.scale = FALSE, tol = 1e-09,  maxit = 500) {
    curve.peeling <- function(x, y) {
        n <- length(y)
        res <- NA
        Fmin <- Inf
        for (i in (n - 3):3) {
            parms1 <- tryCatch(lm(log(y[(i + 1):n]) ~ x[(i + 
                1):n])$coef, error = function(e) rep(NA, 2))
            if (!any(is.na(parms1))) {
                b2 <- parms1[2] * (-1)
                a2 <- exp(parms1[1])
                ynew <- abs(y - a2 * exp(-b2 * x))
                parms2 <- tryCatch(lm(log(ynew[1:i]) ~ x[1:i])$coef, 
                  error = function(e) rep(NA, 2))
                if (!any(is.na(parms2))) {
                  b1 <- parms2[2] * (-1)
                  a1 <- exp(parms2[1])
                  F <- sum((y - (a1 * exp(-b1 * x) + a2 * exp(-b2 * 
                    x))) * (y - (a1 * exp(-b1 * x) + a2 * exp(-b2 * 
                    x))))
                  if (F < Fmin && all(b1 > 0, b2 > 0, b1 > b2)) {
                    res <- as.double(c(a1 = a1, dl = log(b1) - 
                      log(b2), a2 = a2, b2 = log(b2)))
                    Fmin <- F
                  }
                }
            }
        }
        if (any(is.na(res))) {
            parms <- tryCatch(lm(log(y) ~ x), error = function(e) rep(NA, 
                2))
            if (!any(is.na(parms))) {
                b <- parms$coef[2] * (-1)
                a <- exp(parms$coef[1])
                F <- sum(parms$resid * parms$resid)
                if (b > 0) {
                  res <- as.double(c(a = a, b = log(b)))
                  Fmin <- F
                }
            }
        }
        return(res)
    }
    if (!is.vector(time) || !is.vector(conc)) {
        stop("argument time and/or conc invalid")
    }
    if (length(time) != length(conc)) {
        stop("time and conc differ in length")
    }
    if (any(time < 0)) {
        stop("at least one timepoint below zero")
    }
    data <- na.omit(data.frame(conc, time))
    if (any(data$conc <= 0)) {
        data$conc[data$conc <= 0] <- NA
        warning("concentration below or equal to zero were omitted")
        data <- na.omit(data)
    }
    if (nrow(data) < 4) {
        stop("a minimum of 4 observations are required")
    }
    fun <- function(x) {
        return(x)
    }
    if (log.scale) {
        fun <- log
    }
    n <- nrow(data)
    time <- data$time
    conc <- data$conc
    biexploss <- function(par, fun) {
        sum((fun(conc) - fun(par[1] * exp(-(exp(par[4]) + exp(par[2])) * 
            time) + par[3] * exp(-exp(par[4]) * time)))^2)
    }
    singleloss <- function(par, fun) {
        sum((fun(conc) - fun(par[1] * exp(-exp(par[2]) * time)))^2)
    }
    start <- curve.peeling(y = conc, x = time)
    type <- as.character(length(start))
    switch(type, `4` = {
        sum.resid <- biexploss(par = start, fun = fun)
    }, `2` = {
        sum.resid <- singleloss(par = start, fun = fun)
    }, )
    if (sum.resid == Inf) {
        return(warning("no model evaluated"))
    }
    switch(type, `4` = {
        sol <- optim(par = start, fn = biexploss, fun = fun, 
            method = c("Nelder-Mead"), control = list(reltol = tol, 
                maxit = maxit))$par
        b1 <- (exp(sol[4]) + exp(sol[2]))
        a1 <- sol[1]
        b2 <- exp(sol[4])
        a2 <- sol[3]
    }, `2` = {
        sol <- optim(par = start, fn = singleloss, fun = fun, 
            method = c("Nelder-Mead"), control = list(reltol = tol, 
                maxit = maxit))$par
        b1 <- exp(sol[2])
        a1 <- sol[1]
        b2 <- b1
        a2 <- a1
    }, )
    init.hl <- log(2)/b1
    term.hl <- log(2)/b2
    parms <- data.frame(initial = as.double(c(init.hl, b1, a1)), 
        terminal = as.double(c(term.hl, b2, a2)))
    rownames(parms) <- c("halflife", "slope", "intercept")
    res <- list(parms = parms, conc = conc, time = time, method = "biexp")
    class(res) <- "halflife"
    return(res)
}
