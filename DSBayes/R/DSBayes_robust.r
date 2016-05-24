DSBayes <- function (obj = NULL, thetahat, C, lvector, control = list(),
    ...)
{
    ctrl <- list(tol = 0.01, epsilon = 0.005, ci = 0.95, k = NULL,
        transform = "logit", print=TRUE)
    namc <- names(control)
    if (!all(namc %in% names(ctrl)))
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control
    tol <- ctrl$tol
    epsilon <- ctrl$epsilon
    print <- ctrl$print
    ci <- ctrl$ci
    k <- if (is.null(ctrl$k))
        qnorm((9+ci)/10)
    else ctrl$k
    transform <- ctrl$transform
    if (is.null(transform))
        transform <- "none"
    if (print) print(ctrl)
    if (!is.null(obj)) {
        if (names(coef(obj))[1] == "(Intercept)") {
            thetahat <- coef(obj)[-1]
            C <- vcov(obj)[-1, -1]
        }
        else {
            thetahat <- coef(obj)
            C <- vcov(obj)
        }
    }
    p <- length(thetahat)
    m <- (p - 1)/2
    C22 <- C[(m + 2):p, (m + 2):p, drop = FALSE]
    Ceig <- eigen(C22, symmetric = TRUE)
    Q <- Ceig$vec
    Z <- Ceig$value
    Y <- c(crossprod(Q, thetahat[(m + 2):p]))
    Cinv <- solve(C)
    b <- c(Cinv %*% thetahat)
    Binv <- Cinv
    MLE <- function(lvector) {
        nr <- nrow(lvector)
        reg.out <- matrix(NA, nr, 3)
        kf <- qnorm(ci/2 + 0.5)
        for (i in 1:nr) {
            lvec <- lvector[i, ]
            theta0 <- sum(lvec * thetahat)
            sigma <- (t(lvec) %*% C %*% lvec)^0.5
            reg.out[i, 1] <- round(theta0, 3)
            reg.out[i, 2] <- round(theta0 - kf * sigma, 3)
            reg.out[i, 3] <- round(theta0 + kf * sigma, 3)
        }
        colnames(reg.out) <- c("reg.est", "lower CI", "upper CI")
        return(reg.out)
    }
    qfunc.denom <- function(xsi) {
        sapply(xsi, function(xsi) {
            logg <- -m/2 * log(2 * pi) - 1/2 * sum(log(Z + xsi)) -
                1/2 * sum(Y^2/(Z + xsi))
            integrand <- exp(logg)/max(xsi, epsilon)
            integrand
        })
    }
    denom <- integrate(qfunc.denom, lower = 0, upper = Inf)$val
    qfunc.numer <- function(xsi, theta) {
        sapply(xsi, function(xsi) {
            dxsi <- c(rep(0, m + 1), rep(1/xsi, m))
            diag(Binv) <- diag(Cinv) + dxsi
            B <- solve(Binv)
            mean <- t(lvec) %*% B %*% b
            var <- t(lvec) %*% B %*% lvec
            logq <- dnorm(theta, mean = mean, sd = sqrt(var),
                log = TRUE)
            logg <- -m/2 * log(2 * pi) - 1/2 * sum(log(Z + xsi)) -
                1/2 * sum(Y^2/(Z + xsi))
            integrand <- exp(logg + logq)/max(xsi, epsilon)
            integrand
        })
    }
    ds.numer <- function(theta0, thetahat, C, lvec) {
        num <- integrate(qfunc.numer, lower = 0, upper = Inf,
            theta = theta0, rel.tol = tol)$val
        return(num)
    }
    ds.density <- function(theta0, thetahat, C, lvec, denom) sapply(theta0,
        function(x) ds.numer(x, thetahat, C, lvec)/denom)
    klo <- (1 - ci)/2
    khi <- (1 + ci)/2
    if (transform == "logit") {
        lower <- function(x) qlogis(integrate(ds.density, lower = -Inf,
            upper = x, lvec = lvec, thetahat = thetahat, C = C,
            denom = denom, rel.tol = tol)$val) - qlogis(klo)
        upper <- function(x) qlogis(khi) - qlogis(integrate(ds.density,
            lower = -Inf, upper = x, lvec = lvec, thetahat = thetahat,
            C = C, denom = denom, rel.tol = tol)$val)
    }
    else {
        lower <- function(x) integrate(ds.density, lower = -Inf,
            upper = x, lvec = lvec, thetahat = thetahat, C = C,
            denom = denom, rel.tol = tol)$val - klo
        upper <- function(x) khi - integrate(ds.density, lower = -Inf,
            upper = x, lvec = lvec, thetahat = thetahat, C = C,
            denom = denom, rel.tol = tol)$val
    }
    lvector <- as.matrix(lvector)
    if (ncol(lvector) == 1) {
        lvector <- t(lvector)
    }
    nrow <- nrow(lvector)
    out <- matrix(NA, nrow, 3)
    for (i in 1:nrow) {
        if (print) cat(i, "out of ", nrow, "\n")
        lvec <- lvector[i, ]
        theta0 <- sum(lvec * thetahat)
        sigma <- c(t(lvec) %*% C %*% lvec)^0.5
        int.up <- c(theta0 + k * sigma)
        int.low <- c(theta0 - k * sigma)
        mode <- optimize(f = ds.numer, interval = c(int.low,
            int.up), lvec = lvec, thetahat = thetahat, C = C,
            maximum = TRUE, tol = tol)$max
        mod.int.up <- c(mode + k * sigma)
        mod.int.low <- c(mode - k * sigma)

        ci.lower <- try(uniroot(f = lower, interval = c(mod.int.low,
            mode), tol = tol)$root,silent=TRUE)
        if(inherits(ci.lower, "try-error")) ci.lower <- try(dfsane(fn=lower, par=mode,control=list(tol=tol,trace=FALSE))$par,silent=TRUE)
        if(inherits(ci.lower, "try-error")) ci.lower <- NA
        ci.upper <- try(uniroot(f = upper, interval = c(mode, mod.int.up),
            tol = tol)$root,silent=TRUE)
        if(inherits(ci.upper, "try-error")) ci.upper <- try(dfsane(fn=upper, par=mode,control=list(tol=tol,trace=FALSE))$par,silent=TRUE)
        if(inherits(ci.upper, "try-error")) ci.upper <- NA
        out[i, ] <- round(c(mode, ci.lower, ci.upper), 3)
    }
    colnames(out) <- c("Mode", "lower CI", "upper CI")
    mle <- MLE(lvector)
    list(bayes = out, mle = mle)
}
