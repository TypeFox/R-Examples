SN1<- function (mu.link = "identity", sigma.link = "log", nu.link = "identity") 
{
    mstats <- checklink("mu.link", "Skew normal type 1 (Azzalini type 1)", 
        substitute(mu.link), c("inverse", "log", "identity", 
            "own"))
    dstats <- checklink("sigma.link", "Skew normal type 1 (Azzalini type 1)", 
        substitute(sigma.link), c("inverse", "log", "identity", 
            "own"))
    vstats <- checklink("nu.link", "Skew normal type 1 (Azzalini type 1)", 
        substitute(nu.link), c("inverse", "log", "identity", 
            "own"))
    structure(list(family = c("SN1", "Skew normal type 1 (Azzalini type 1)"), 
        parameters = list(mu = TRUE, sigma = TRUE, nu = TRUE), nopar = 3, 
                   type = "Continuous", mu.link = as.character(substitute(mu.link)), 
        sigma.link = as.character(substitute(sigma.link)), nu.link = as.character(substitute(nu.link)), 
        mu.linkfun = mstats$linkfun, 
        sigma.linkfun = dstats$linkfun, nu.linkfun = vstats$linkfun, 
        mu.linkinv = mstats$linkinv, 
        sigma.linkinv = dstats$linkinv, nu.linkinv = vstats$linkinv, 
        mu.dr = mstats$mu.eta, 
        sigma.dr = dstats$mu.eta, nu.dr = vstats$mu.eta, 
        dldm = function(y, mu, sigma, nu) {
            z <- (y - mu)/sigma
            w <- nu * z
            s <- ((abs(w))^2)/2
            lpdf <- (1 - (1/2)) * log(2) - s - lgamma(1/2) - 
                log(2)
            lcdf <- log(0.5 * (1 + pgamma(s, shape = 1/2, scale = 1) * 
                sign(w)))
            dldm <- -(exp(lpdf - lcdf)) * nu/sigma + sign(z) * 
                (abs(z)^(2 - 1))/sigma
            dldm
        }, d2ldm2 = function(y, mu, sigma, nu) {
            z <- (y - mu)/sigma
            w <- nu * z
            s <- ((abs(w))^2)/2
            lpdf <- (1 - (1/2)) * log(2) - s - lgamma(1/2) - 
                log(2)
            lcdf <- log(0.5 * (1 + pgamma(s, shape = 1/2, scale = 1) * 
                sign(w)))
            dldm <- -(exp(lpdf - lcdf)) * nu/sigma + sign(z) * 
                (abs(z)^(2 - 1))/sigma
            d2ldm2 <- -dldm * dldm
            d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2, -1e-15)
            d2ldm2
        }, dldd = function(y, mu, sigma, nu) {
            z <- (y - mu)/sigma
            w <- nu * z
            s <- ((abs(w))^2)/2
            lpdf <- (1 - (1/2)) * log(2) - s - lgamma(1/2) - 
                log(2)
            lcdf <- log(0.5 * (1 + pgamma(s, shape = 1/2, scale = 1) * 
                sign(w)))
            dldd <- -(exp(lpdf - lcdf)) * nu * z/sigma + ((abs(z)^2) - 
                1)/sigma
            dldd
        }, d2ldd2 = function(y, mu, sigma, nu) {
            z <- (y - mu)/sigma
            w <- nu * z
            s <- ((abs(w))^2)/2
            lpdf <- (1 - (1/2)) * log(2) - s - lgamma(1/2) - 
                log(2)
            lcdf <- log(0.5 * (1 + pgamma(s, shape = 1/2, scale = 1) * 
                sign(w)))
            dldd <- -(exp(lpdf - lcdf)) * nu * z/sigma + ((abs(z)^2) - 
                1)/sigma
            d2ldd2 <- -dldd * dldd
            d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2, -1e-15)
            d2ldd2
        }, dldv = function(y, mu, sigma, nu) {
            z <- (y - mu)/sigma
            w <- nu * z
            dwdv <- w/nu
            s <- ((abs(w))^2)/2
            lpdf <- (1 - (1/2)) * log(2) - s - lgamma(1/2) - 
                log(2)
            lcdf <- log(0.5 * (1 + pgamma(s, shape = 1/2, scale = 1) * 
                sign(w)))
            dldv <- (exp(lpdf - lcdf)) * dwdv
            dldv
        }, d2ldv2 = function(y, mu, sigma, nu) {
            z <- (y - mu)/sigma
            w <- nu * z
            dwdv <- w/nu
            s <- ((abs(w))^2)/2
            lpdf <- (1 - (1/2)) * log(2) - s - lgamma(1/2) - 
                log(2)
            lcdf <- log(0.5 * (1 + pgamma(s, shape = 1/2, scale = 1) * 
                sign(w)))
            dldv <- (exp(lpdf - lcdf)) * dwdv
            d2ldv2 <- -dldv * dldv
            d2ldv2 <- ifelse(d2ldv2 < -1e-15, d2ldv2, -1e-15)
            d2ldv2
        }, d2ldmdd = function(y, mu, sigma, nu) {
            z <- (y - mu)/sigma
            w <- nu * z
            s <- ((abs(w))^2)/2
            lpdf <- (1 - (1/2)) * log(2) - s - lgamma(1/2) - 
                log(2)
            lcdf <- log(0.5 * (1 + pgamma(s, shape = 1/2, scale = 1) * 
                sign(w)))
            dldm <- -(exp(lpdf - lcdf)) * nu/sigma + sign(z) * 
                (abs(z)^(2 - 1))/sigma
            dldd <- -(exp(lpdf - lcdf)) * nu * z/sigma + ((abs(z)^2) - 
                1)/sigma
            d2ldmdd <- -(dldm * dldd)
            d2ldmdd
        }, d2ldmdv = function(y, mu, sigma, nu) {
            z <- (y - mu)/sigma
            w <- nu * z
            s <- ((abs(w))^2)/2
            lpdf <- (1 - (1/2)) * log(2) - s - lgamma(1/2) - 
                log(2)
            lcdf <- log(0.5 * (1 + pgamma(s, shape = 1/2, scale = 1) * 
                sign(w)))
            dldm <- -(exp(lpdf - lcdf)) * nu/sigma + sign(z) * 
                (abs(z)^(2 - 1))/sigma
            dwdv <- w/nu
            dldv <- (exp(lpdf - lcdf)) * dwdv
            d2ldmdv <- -(dldm * dldv)
            d2ldmdv
        }, d2ldddv = function(y, mu, sigma, nu) {
            z <- (y - mu)/sigma
            w <- nu * z
            s <- ((abs(w))^2)/2
            lpdf <- (1 - (1/2)) * log(2) - s - lgamma(1/2) - 
                log(2)
            lcdf <- log(0.5 * (1 + pgamma(s, shape = 1/2, scale = 1) * 
                sign(w)))
            dldd <- -(exp(lpdf - lcdf)) * nu * z/sigma + ((abs(z)^2) - 
                1)/sigma
            dwdv <- w/nu
            dldv <- (exp(lpdf - lcdf)) * dwdv
            d2ldddv <- -(dldd * dldv)
            d2ldddv
        }, G.dev.incr = function(y, mu, sigma, nu, ...) {
            -2 * dSN1(y, mu, sigma, nu, log = TRUE)
        }, rqres = expression(rqres(pfun = "pSN1", type = "Continuous", 
            y = y, mu = mu, sigma = sigma, nu = nu)), 
        mu.initial = expression(mu <- (y + mean(y))/2), sigma.initial = expression(sigma <- rep(sd(y)/4, 
            length(y))), nu.initial = expression(nu <- rep(0.1, 
            length(y))), mu.valid = function(mu) TRUE, sigma.valid = function(sigma) all(sigma > 
            0), nu.valid = function(nu) TRUE, y.valid = function(y) TRUE), class = c("gamlss.family", 
        "family"))
}
                   
#dSN1
dSN1<-function (x, mu = 0, sigma = 1, nu = 0, log = FALSE) 
{
    if (any(sigma < 0)) 
        stop(paste("sigma must be positive", "\n", ""))
    z <- (x - mu)/sigma
    w <- nu * z
    sz <- ((abs(z))^2)/2
    s <- ((abs(w))^2)/2
    lpdf <- (1 - (1/2)) * log(2) - sz - lgamma(1/2) - log(2)
    lcdf1 <- log(0.5 * (1 + pgamma(s, shape = 1/2, scale = 1) * 
        sign(w)))
    cdf2 <- 0.5 + w * exp((1 - (1/2)) * log(2) - lgamma(1/2) - 
        log(2))
    suppressWarnings(lcdf2 <- log(cdf2))
    lcdf <- ifelse((s == 0), lcdf2, lcdf1)
    loglik <- lpdf + lcdf + log(2) - log(sigma)
    if (log == FALSE) 
        ft <- exp(loglik)
    else ft <- loglik
    ft
}

#pSN1
pSN1<-function (q, mu = 0, sigma = 1, nu = 0, lower.tail = TRUE, 
    log.p = FALSE) 
{
    if (any(sigma <= 0)) 
        stop(paste("sigma must be positive", "\n", ""))
    lp <- pmax.int(length(q), length(mu), length(sigma), length(nu))
    q <- rep(q, length = lp)
    sigma <- rep(sigma, length = lp)
    mu <- rep(mu, length = lp)
    nu <- rep(nu, length = lp)
    cdf <- rep(0, length = lp)
    for (i in 1:lp) {
        cdf[i] <- integrate(function(x) dSN1(x, mu = 0, sigma = 1, 
            nu = nu[i]), -Inf, (q[i] - mu[i])/sigma[i])$value
    }
    if (lower.tail == TRUE) 
        cdf <- cdf
    else cdf <- 1 - cdf
    if (log.p == FALSE) 
        cdf <- cdf
    else cdf <- log(cdf)
    cdf
}
      
#qSN1
qSN1<-function (p, mu = 0, sigma = 1, nu = 0, lower.tail = TRUE, 
    log.p = FALSE) 
{
    h1 <- function(q) {
        pSN1(q, mu = mu[i], sigma = sigma[i], nu = nu[i]) - 
            p[i]
    }
    h <- function(q) {
        pSN1(q, mu = mu[i], sigma = sigma[i], nu = nu[i])
    }
    if (any(sigma <= 0)) 
        stop(paste("sigma must be positive", "\n", ""))
    if (log.p == TRUE) 
        p <- exp(p)
    else p <- p
    if (lower.tail == TRUE) 
        p <- p
    else p <- 1 - p
    if (any(p < 0) | any(p > 1)) 
        stop(paste("p must be between 0 and 1", "\n", ""))
    lp <- max(length(p), length(mu), length(sigma), length(nu))
    p <- rep(p, length = lp)
    sigma <- rep(sigma, length = lp)
    mu <- rep(mu, length = lp)
    nu <- rep(nu, length = lp)
    q <- rep(0, lp)
    for (i in seq(along = p)) {
        if (h(mu[i]) < p[i]) {
            interval <- c(mu[i], mu[i] + sigma[i])
            j <- 2
            while (h(interval[2]) < p[i]) {
                interval[2] <- mu[i] + j * sigma[i]
                j <- j + 1
            }
        }
        else {
            interval <- c(mu[i] - sigma[i], mu[i])
            j <- 2
            while (h(interval[1]) > p[i]) {
                interval[1] <- mu[i] - j * sigma[i]
                j <- j + 1
            }
        }
        q[i] <- uniroot(h1, interval)$root
    }
    q
}
        
rSN1<-function(n, mu = 0, sigma = 1, nu = 0) 
{
    if (any(sigma <= 0)) 
        stop(paste("sigma must be positive", "\n", ""))
    n <- ceiling(n)
    p <- runif(n)
    r <- qSN1(p, mu = mu, sigma = sigma, nu = nu)
    r
}