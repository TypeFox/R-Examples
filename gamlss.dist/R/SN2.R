###Skew Normal Type 2

#SN2
SN2<-function (mu.link = "identity", sigma.link = "log", nu.link = "log") 
{
    mstats <- checklink("mu.link", "skew normal type 2", 
        substitute(mu.link), c("inverse", "log", "identity", 
            "own"))
    dstats <- checklink("sigma.link", "skew normal type 2", 
        substitute(sigma.link), c("inverse", "log", "identity", 
            "own"))
    vstats <- checklink("nu.link", "skew normal type 2", 
        substitute(nu.link), c("inverse", "log", "identity", 
            "own"))
    structure(list(family = c("SN2", "skew normal type 2"), 
        parameters = list(mu = TRUE, sigma = TRUE, nu = TRUE), nopar = 3, type = "Continuous", mu.link = as.character(substitute(mu.link)), 
        sigma.link = as.character(substitute(sigma.link)), nu.link = as.character(substitute(nu.link)), 
        mu.linkfun = mstats$linkfun, 
        sigma.linkfun = dstats$linkfun, nu.linkfun = vstats$linkfun, 
        mu.linkinv = mstats$linkinv, 
        sigma.linkinv = dstats$linkinv, nu.linkinv = vstats$linkinv, 
        mu.dr = mstats$mu.eta, 
        sigma.dr = dstats$mu.eta, nu.dr = vstats$mu.eta, 
        dldm = function(y, mu, sigma, nu) {
            z <- (y - mu)/sigma
            dldma <- sign(z) * (nu * 2/(2 * sigma)) * ((nu * 
                abs(z))^(2 - 1))
            dldmb <- sign(z) * (2/(2 * sigma * nu)) * ((abs(z)/nu)^(2 - 
                1))
            dldm <- ifelse(y < mu, dldma, dldmb)
            dldm
        }, d2ldm2 = function(y, mu, sigma, nu) {
            z <- (y - mu)/sigma
            dldma <- sign(z) * (nu * 2/(2 * sigma)) * ((nu * 
                abs(z))^(2 - 1))
            dldmb <- sign(z) * (2/(2 * sigma * nu)) * ((abs(z)/nu)^(2 - 
                1))
            dldm <- ifelse(y < mu, dldma, dldmb)
            d2ldm2 <- -dldm * dldm
            d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2, -1e-15)
            d2ldm2
        }, dldd = function(y, mu, sigma, nu) {
            z <- (y - mu)/sigma
            dldda <- (nu * abs(z))^2
            dlddb <- (abs(z)/nu)^2
            dldd <- ifelse(y < mu, dldda, dlddb)
            dldd <- dldd * (2/(2 * sigma)) - (1/sigma)
            dldd
        }, d2ldd2 = function(y, mu, sigma, nu) {
            z <- (y - mu)/sigma
            dldda <- (nu * abs(z))^2
            dlddb <- (abs(z)/nu)^2
            dldd <- ifelse(y < mu, dldda, dlddb)
            dldd <- dldd * (2/(2 * sigma)) - (1/sigma)
            d2ldd2 <- -dldd * dldd
            d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2, -1e-15)
            d2ldd2
        }, dldv = function(y, mu, sigma, nu) {
            z <- (y - mu)/sigma
            dldva <- sign(z) * ((nu * abs(z))^2)
            dldvb <- sign(z) * ((abs(z)/nu)^2)
            dldv <- ifelse(y < mu, dldva, dldvb)
            dldv <- dldv * (2/(2 * nu)) + (1/nu) - (2 * nu)/(1 + 
                (nu^2))
            dldv
        }, d2ldv2 = function(y, mu, sigma, nu) {
            z <- (y - mu)/sigma
            dldva <- sign(z) * ((nu * abs(z))^2)
            dldvb <- sign(z) * ((abs(z)/nu)^2)
            dldv <- ifelse(y < mu, dldva, dldvb)
            dldv <- dldv * (2/(2 * nu)) + (1/nu) - (2 * nu)/(1 + 
                (nu^2))
            d2ldv2 <- -dldv * dldv
            d2ldv2 <- ifelse(d2ldv2 < -1e-04, d2ldv2, -1e-04)
            d2ldv2
        }, d2ldmdd = function(y, mu, sigma, nu) {
            z <- (y - mu)/sigma
            dldma <- sign(z) * (nu * 2/(2 * sigma)) * ((nu * 
                abs(z))^(2 - 1))
            dldmb <- sign(z) * (2/(2 * sigma * nu)) * ((abs(z)/nu)^(2 - 
                1))
            dldm <- ifelse(y < mu, dldma, dldmb)
            dldda <- (nu * abs(z))^2
            dlddb <- (abs(z)/nu)^2
            dldd <- ifelse(y < mu, dldda, dlddb)
            dldd <- dldd * (2/(2 * sigma)) - (1/sigma)
            d2ldmdd <- -(dldm * dldd)
            d2ldmdd
        }, d2ldmdv = function(y, mu, sigma, nu) {
            z <- (y - mu)/sigma
            dldma <- sign(z) * (nu * 2/(2 * sigma)) * ((nu * 
                abs(z))^(2 - 1))
            dldmb <- sign(z) * (2/(2 * sigma * nu)) * ((abs(z)/nu)^(2 - 
                1))
            dldm <- ifelse(y < mu, dldma, dldmb)
            dldva <- sign(z) * ((nu * abs(z))^2)
            dldvb <- sign(z) * ((abs(z)/nu)^2)
            dldv <- ifelse(y < mu, dldva, dldvb)
            dldv <- dldv * (2/(2 * nu)) + (1/nu) - (2 * nu)/(1 + 
                (nu^2))
            d2ldmdv <- -(dldm * dldv)
            d2ldmdv
        }, d2ldddv = function(y, mu, sigma, nu) {
            z <- (y - mu)/sigma
            dldda <- (nu * abs(z))^2
            dlddb <- (abs(z)/nu)^2
            dldd <- ifelse(y < mu, dldda, dlddb)
            dldd <- dldd * (2/(2 * sigma)) - (1/sigma)
            dldva <- sign(z) * ((nu * abs(z))^2)
            dldvb <- sign(z) * ((abs(z)/nu)^2)
            dldv <- ifelse(y < mu, dldva, dldvb)
            dldv <- dldv * (2/(2 * nu)) + (1/nu) - (2 * nu)/(1 + 
                (nu^2))
            d2ldddv <- -(dldd * dldv)
            d2ldddv
        }, G.dev.incr = function(y, mu, sigma, nu, ...) -2 * 
            dSN2(y, mu, sigma, nu, log = TRUE), rqres = expression(rqres(pfun = "pSN2", 
            type = "Continuous", y = y, mu = mu, sigma = sigma, 
            nu = nu)), mu.initial = expression(mu <- (y + 
            mean(y))/2), sigma.initial = expression(sigma <- rep(sd(y), 
            length(y))), nu.initial = expression(nu <- rep(1, 
            length(y))), mu.valid = function(mu) TRUE, sigma.valid = function(sigma) all(sigma > 
            0), nu.valid = function(nu) all(nu > 0), y.valid = function(y) TRUE), 
            class = c("gamlss.family", "family"))
}

#dSN2
dSN2<-function (x, mu = 0, sigma = 1, nu = 2, log = FALSE)
{
    if (any(sigma <= 0)) 
        stop(paste("sigma must be positive", "\n", ""))
    if (any(nu <= 0)) 
        stop(paste("nu must be positive", "\n", ""))
    z <- (x - mu)/sigma
    suppressWarnings(loglik1 <- -0.5 * ((nu * abs(z))^2))
    suppressWarnings(loglik2 <- -0.5 * ((abs(z)/nu)^2))
    loglik <- ifelse(x < mu, loglik1, loglik2)
    loglik <- loglik - log(sigma) + log(nu) - log(1 + (nu^2)) - 
        (1/2) * log(2) - lgamma(1 + (1/2))
    fy <- if (log == FALSE) 
        exp(loglik)
    else loglik
    fy
}
    
#pSN2
pSN2<-function (q, mu = 0, sigma = 1, nu = 2, lower.tail = TRUE, 
    log.p = FALSE) 
{
    if (any(sigma < 0)) 
        stop(paste("sigma must be positive", "\n", ""))
    if (any(nu <= 0)) 
        stop(paste("nu must be positive", "\n", ""))
    k <- nu^2
    z1 <- nu * (q - mu)/(sigma * (2^(1/2)))
    z2 <- (q - mu)/(sigma * nu * (2^(1/2)))
    s1 <- (abs(z1)^2)
    s2 <- (abs(z2)^2)
    cdf1 <- 1 - pgamma(s1, shape = 1/2, scale = 1)
    cdf2 <- 1 + k * pgamma(s2, shape = 1/2, scale = 1)
    cdf <- ifelse(q < mu, cdf1, cdf2)
    cdf <- cdf/(1 + k)
    if (length(2) > 1) 
        cdf <- ifelse(2 > 10000, (q - (mu - (sigma/nu)))/(sigma * 
            ((1/nu) + nu)), cdf)
    else cdf <- if (2 > 10000) 
        (q - (mu - (sigma/nu)))/(sigma * ((1/nu) + nu))
    else cdf
    if (lower.tail == TRUE) 
        cdf <- cdf
    else cdf <- 1 - cdf
    if (log.p == FALSE) 
        cdf <- cdf
    else cdf <- log(cdf)
    cdf
}
    
#qSN2
qSN2<-function (p, mu = 0, sigma = 1, nu = 2, lower.tail = TRUE, 
    log.p = FALSE) 
{
    if (any(sigma < 0)) 
        stop(paste("sigma must be positive", "\n", ""))
    if (any(nu <= 0)) 
        stop(paste("nu must be positive", "\n", ""))
    if (log.p == TRUE) 
        p <- exp(p)
    else p <- p
    if (any(p <= 0) | any(p >= 1)) 
        stop(paste("p must be between 0 and 1", "\n", ""))
    if (lower.tail == TRUE) 
        p <- p
    else p <- 1 - p
    k <- nu^2
    suppressWarnings(q1 <- mu - (sigma * (2^(1/2))/nu) * ((qgamma(1 - 
        p * (1 + k), shape = 1/2, scale = 1))^(1/2)))
    suppressWarnings(q2 <- mu + (sigma * nu * (2^(1/2))) * 
        ((qgamma((-1/k) * (1 - p * (1 + k)), shape = 1/2, scale = 1))^(1/2)))
    q <- ifelse(p < (1/(1 + k)), q1, q2)
    q
}

#rSN2
rSN2<-function (n, mu = 0, sigma = 1, nu = 2) 
{
    if (any(sigma <= 0)) 
        stop(paste("sigma must be positive", "\n", ""))
    if (any(nu <= 0)) 
        stop(paste("nu must be positive", "\n", ""))
    if (any(n <= 0)) 
        stop(paste("n must be a positive integer", "\n", ""))
    n <- ceiling(n)
    p <- runif(n)
    r <- qSN2(p, mu = mu, sigma = sigma, nu = nu)
    r
}
