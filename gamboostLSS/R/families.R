###
# Constructor Function

Families <- function(..., qfun = NULL, name = NULL) {
    RET <- list(...)
    class(RET) <- "families"
    ## check if response function is not specified
    check <- sapply(RET, function(x){
        bdy <- body(x@response)
        (length(x@response(c(-1,0,1))) != 3 && class(bdy) != "call" &&
         length(bdy) == 1 && is.na(bdy))
    })
    if (any(check))
        stop("response function not specified in families for:\n\t",
             paste(names(RET)[check], collapse =", "))
    attr(RET, "qfun") <- qfun
    attr(RET, "name") <- name
    return(RET)
}

###
# Negative Binomial LSS Family

NBinomialMu <- function(mu = NULL, sigma = NULL, stabilization) {
    loss <- function(sigma, y, f)
        -(lgamma(y + sigma) - lgamma(sigma) - lgamma(y + 1) + sigma * log(sigma)
          + y * f - (y + sigma) * log(exp(f) + sigma))
    risk <- function(y, f, w = 1){
        RET <- sum(w * loss(y = y, f = f, sigma = sigma))
        return(RET)
    }
    ngradient <- function(y, f, w = 1) {
        ngr <- y - (y + sigma)/(exp(f) + sigma) * exp(f)
        ngr <- stabilize_ngradient(ngr, w = w, stabilization)
        return(ngr)
    }
    offset <- function(y, w){
        if (!is.null(mu)){
            RET <- log(mu)
        } else {
            if (is.null(sigma)) {
                sigma <<- mean(y)^2 / (var(y) - mean(y))
                sigma <<- ifelse(sigma < 0, 1e-10, sigma)
            }
            ### look for starting value of f = log(mu) in "interval"
            ### i.e. mu possibly ranges from 1e-10 to 1e10
            RET <- optimize(risk, interval = c(log(1e-10), log(1e10)), y = y, w = w)$minimum
        }
        return(RET)
    }

    Family(ngradient = ngradient, offset = offset, risk = risk, loss = loss,
           response = function(f) exp(f),
           check_y = function(y) {
               if (!all.equal(unique(y - floor(y)), 0))
                   stop("response must be an integer")
               if (any(y < 0))
                   stop("response must be >= 0")
               y
           }, name = "Negative Negative-Binomial Likelihood: mu (log link)")
}

NBinomialSigma <- function(mu = NULL, sigma = NULL, stabilization) {
    # this family boosts sigma therefore f is sigma
    loss <- function(mu, y, f)
        -(lgamma(y + exp(f)) - lgamma(exp(f)) - lgamma(y + 1) + exp(f) * f +
          y * log(mu) - (y + exp(f)) * log(mu + exp(f)))
    risk <- function(y, f, w = 1){
        RET <- sum(w * loss(y = y, f = f, mu = mu))
        return(RET)
    }
    ngradient <- function(y, f, w = 1) {       # f is sigma !
        ngr <- exp(f)*(digamma(y +exp(f)) - digamma(exp(f)) + log(exp(f)) + 1 -
                       log(mu +exp(f)) - (exp(f) + y)/(mu +exp(f)))
        ngr <- stabilize_ngradient(ngr, w = w, stabilization)
        return(ngr)
    }
    offset <- function(y,w){
        if (!is.null(sigma)){
            RET <- log(sigma)
        } else {
            if (is.null(mu))
                mu <<- mean(y)
            ### look for starting value of f = log(sigma) in "interval"
            ### i.e. sigma possibly ranges from 1e-10 to 1e10
            RET <- optimize(risk, interval = c(log(1e-10), log(1e10)), y = y, w = w)$minimum
        }
        return(RET)
    }

    Family(ngradient = ngradient, offset = offset, risk = risk, loss = loss,
           response = function(f) exp(f),
           check_y = function(y) {
               if (!all.equal(unique(y - floor(y)), 0))
                   stop("response must be an integer")
               if (any(y < 0))
                   stop("response must be >= 0")
               y
           }, name = "Negative Negative-Binomial Likelihood: sigma (log link)")
}


NBinomialLSS <- function(mu = NULL, sigma = NULL,
                         stabilization = c("none", "MAD")) {
    if ((!is.null(sigma) && sigma <= 0) || (!is.null(mu) && mu <= 0))
        stop(sQuote("sigma"), " and ", sQuote("mu"),
             " must be greater than zero")
    stabilization <- check_stabilization(stabilization)
    Families(mu = NBinomialMu(mu = mu, sigma = sigma, stabilization = stabilization),
             sigma = NBinomialSigma(mu = mu, sigma = sigma, stabilization = stabilization),
             qfun = qNBinomial,
             name = "Negative Binomial")
}

## we use almost the same parameterization as NBI but with theta = 1/sigma
qNBinomial <- function(p, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
    ## require gamlss.dist
    if (!requireNamespace("gamlss.dist", quietly = TRUE))
        stop("Please install package 'gamlss.dist' for using qNBinomial.")
    gamlss.dist::qNBI(p = p, mu = mu, sigma = 1/sigma, lower.tail = lower.tail, log.p = log.p)
}


###
# T-distribution LSS Family

StudentTMu <- function(mu = NULL, sigma = NULL, df = NULL, stabilization) {
    loss <- function(df, sigma,y, f){
        -1 * (lgamma((df+1)/2) - log(sigma) - lgamma(1/2) - lgamma(df/2) - 0.5 *
              log(df) - (df+1)/2 * log(1 + (y-f)^2 / (df * sigma^2)))
    }
    risk <- function(y, f, w = 1){
        sum(w * loss(y = y, f = f, df = df, sigma = sigma))
    }
    ngradient <- function(y, f, w = 1) {
        ngr <- (df+1)*(y-f)/(df*sigma^2 +(y-f)^2)
        ngr <- stabilize_ngradient(ngr, w = w, stabilization)
        return(ngr)
     }

    offset <- function(y, w){
        if (!is.null(mu)){
            RET <- mu
        } else {
            if (is.null(sigma))
                sigma <<- 1
            if (is.null(df))
                df <<- 4
            RET <- optimize(risk, interval = range(y), y = y, w = w)$minimum
        }
        return(RET)
    }

    Family(ngradient = ngradient, risk = risk, loss = loss,
           offset=offset,
           response = function(f) f,
           check_y = function(y) {
               if (!is.numeric(y) || !is.null(dim(y)))
                   stop("response is not a numeric vector ")
               y
           },  name = "Student's t-distribution: mu (id link)")
}

StudentTSigma <- function(mu = NULL, sigma = NULL, df = NULL, stabilization) {
    loss <- function(df, mu, y, f){
        -1 * (lgamma((df+1)/2) - f - lgamma(1/2) - lgamma(df/2) - 0.5 * log(df) -
              (df+1)/2 * log(1 + (y-mu)^2 / (df * exp(2 * f))))
    }
    risk <- function(y, f, w = 1){
        sum(w * loss(y = y, f = f, df = df, mu = mu))
    }
    ngradient <- function(y, f, w = 1) {
        ngr <- (-1 + (df+1)/(df*exp(2*f)/(y-mu)^2 + 1))
        ngr <- stabilize_ngradient(ngr, w = w, stabilization)
        return(ngr)
    }
    offset <- function(y, w){
        if (!is.null(sigma)){
            RET <- log(sigma)
        } else {
            if (is.null(mu))
                mu <<- mean(y)
            if (is.null(df))
                df <<- 4
            ### look for starting value of f = log(sigma) in "interval"
            ### i.e. sigma possibly ranges from 1e-10 to 1e10
            RET <- optimize(risk, interval = c(log(1e-10), log(1e10)), y = y, w = w)$minimum
        }
        return(RET)
    }

    Family(ngradient = ngradient, risk = risk, loss = loss,
           offset=offset,
           response = function(f) exp(f),
           check_y = function(y) {
               if (!is.numeric(y) || !is.null(dim(y)))
                   stop("response is not a numeric vector ")
               y
           },  name = "Student's t-distribution: sigma (log link)")
}

StudentTDf <- function(mu = NULL, sigma = NULL, df = NULL, stabilization) {
    loss <- function(sigma, mu,y, f){
        -1 * (lgamma((exp(f)+1)/2) - log(sigma) - lgamma(1/2) - lgamma(exp(f)/2) -
              0.5*f - (exp(f)+1)/2 * log(1 + (y-mu)^2 / (exp(f)*sigma^2)))
    }
    risk <- function(y, f, w = 1){
        sum(w * loss(y = y, f = f, sigma = sigma, mu = mu))
    }
    ngradient <- function(y, f, w = 1) {
        ngr <- exp(f)/2 * (digamma((exp(f)+1)/2)-digamma(exp(f)/2)) - 0.5 -
            (exp(f)/2 * log(1+ (y-mu)^2/(exp(f)*sigma^2)) -
             (y-mu)^2/(sigma^2 + (y-mu)^2/exp(f)) * (exp(-f) +1)/2 )
        ngr <- stabilize_ngradient(ngr, w = w, stabilization)
        return(ngr)
    }
    offset <- function(y, w){
        if (!is.null(df)){
            RET <- log(df)
        } else {
            if (is.null(mu))
                mu <<- mean(y)
            if (is.null(sigma))
                sigma <<- 1
            ### look for starting value of f = log(df) in "interval"
            ### i.e. df possibly ranges from 1e-10 to 1e10
            RET <- optimize(risk, interval = c(log(1e-10), log(1e10)), y = y, w = w)$minimum
        }
        return(RET)
    }

    Family(ngradient = ngradient, risk = risk, loss = loss,
           offset=offset,
           response = function(f) exp(f),
           check_y = function(y) {
               if (!is.numeric(y) || !is.null(dim(y)))
                   stop("response is not a numeric vector ")
               y
           },  name = "Student's t-distribution: df (log link)")
}



StudentTLSS <- function(mu = NULL, sigma = NULL, df = NULL,
                        stabilization = c("none", "MAD")) {
    if ((!is.null(sigma) && sigma <= 0) || (!is.null(df) && df <= 0))
        stop(sQuote("sigma"), " and ", sQuote("df"),
             " must be greater than zero")
    stabilization <- check_stabilization(stabilization)
    Families(mu = StudentTMu(mu = mu, sigma = sigma, df = df, stabilization = stabilization),
             sigma = StudentTSigma(mu = mu, sigma = sigma, df = df, stabilization = stabilization),
             df = StudentTDf(mu = mu, sigma = sigma, df = df, stabilization = stabilization),
             qfun = qT,
             name = "Student T")
}

qT <- function(p, mu, sigma, df, n, lower.tail = TRUE, log.p = FALSE) {
    if (any(sigma <= 0))
        stop("sigma must be greater 0")
    if (any(df <= 0))
        stop("df must be greater 0")
    if (length(n) != 1 || n <= 0)
        stop("n must be a single value greater 0")
    ncp <- mu * sqrt(n)/sigma
    q <- qt(p, df = df, ncp = ncp, lower.tail = lower.tail, log.p = log.p)
    return(q)
}


###
# Log-Normal LSS Family

LogNormalMu <- function (mu = NULL, sigma = NULL, stabilization){
    loss <- function(sigma, y, f) {
        logfw <- function(pred)
            dnorm(pred, log = TRUE)
        logSw <- function(pred)
            pnorm(pred, lower.tail = FALSE, log.p = TRUE)
        eta <- (log(y[,1]) - f)/sigma
        -y[,2] * (logfw(eta) - log(sigma)) - (1 - y[,2]) * logSw(eta)
    }
    risk <- function(y, f, w = 1)
        sum(w * loss(y = y, f = f, sigma = sigma))
    ngradient <- function(y, f, w = 1) {
        eta <- (log(y[,1]) - f)/sigma
        ngr <- (y[,2] * eta + (1 - y[,2]) * dnorm(eta)/(1 - pnorm(eta)))/sigma
        ngr <- stabilize_ngradient(ngr, w = w, stabilization)
        return(ngr)
    }
    offset <- function(y, w){
        if (!is.null(mu)){
            RET <- mu
        } else {
            if (is.null(sigma))
                sigma <<- 1
            RET <- optimize(risk, interval = c(log(1e-10), log(1e10)),
                            y = y, w = w)$minimum
        }
        return(RET)
    }

    Family(ngradient = ngradient, risk = risk, offset = offset, loss = loss,
           response = function(f) f,
           check_y = function(y) {
               if (!inherits(y, "Surv"))
                   stop("response is not an object of class ", sQuote("Surv"),
                        " but ", sQuote("families = LogNormalLSS()"))
               y
           }, name = "Log-Normal AFT Model: mu (id link)")
}

LogNormalSigma <- function(mu = NULL, sigma = NULL, stabilization){
    loss <- function(mu, y, f) {
        logfw <- function(pred)
            dnorm(pred, log = TRUE)
        logSw <- function(pred)
            pnorm(pred, lower.tail = FALSE, log.p = TRUE)
        eta <- (log(y[,1]) - mu) / exp(f)
        -y[,2] * (logfw(eta) - f) - (1 - y[,2]) * logSw(eta)
    }
    risk <- function(y, f, w = 1)
        sum(w * loss(y = y, f = f, mu = mu))
    ngradient <- function(y, f, w = 1) {
        eta <- (log(y[,1]) - mu)/exp(f)
        ngr <- -(y[,2] - y[,2]*eta^2 + (y[,2]-1)*eta*dnorm(eta)/(1-pnorm(eta)))
        ngr <- stabilize_ngradient(ngr, w = w, stabilization)
        return(ngr)
    }
    offset <- function(y, w){
        if (!is.null(sigma)){
            RET <- log(sigma)
        } else {
            if (is.null(mu))
                mu <<- 0
            ## look for starting value of f = log(sigma) in "interval"
            ## i.e. sigma possibly ranges from 1e-10 to 1e10
            RET <- optimize(risk, interval = c(log(1e-10), log(1e10)),
                            y = y, w = w)$minimum
        }
        return(RET)
    }

    Family(ngradient = ngradient, risk = risk, offset = offset, loss = loss,
           response = function(f) exp(f),
           check_y = function(y) {
               if (!inherits(y, "Surv"))
                   stop("response is not an object of class ", sQuote("Surv"),
                        " but ", sQuote("family = LogNormalLSS()"))
               y
           }, name = "Log-Normal AFT Model: sigma (log link)")
}

LogNormalLSS <- function(mu = NULL, sigma = NULL,
                         stabilization = c("none", "MAD")) {
    if ((!is.null(sigma) && sigma <= 0))
        stop(sQuote("sigma"), " must be greater than zero")
    stabilization <- check_stabilization(stabilization)
    Families(mu = LogNormalMu(mu = mu, sigma = sigma, stabilization = stabilization),
             sigma = LogNormalSigma(mu = mu, sigma = sigma, stabilization = stabilization),
             qfun = qLogNormal,
             name = "Log-Normal")
}

qLogNormal <- function(p, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
    qlnorm(p = p, meanlog = log(mu), sdlog = log(sigma), lower.tail = lower.tail,
           log.p = log.p)
}


###
# LogLog LSS Family

LogLogMu <- function (mu = NULL, sigma = NULL, stabilization){
    loss <- function(sigma, y, f) {
        logfw <- function(pred)
            dlogis(pred, log = TRUE)
            #pred - 2 * (1 + exp(pred))
        logSw <- function(pred)
            plogis(pred, lower.tail = FALSE, log.p = TRUE)
            #1/(1 + exp(pred))
        eta <- (log(y[,1]) - f)/sigma
        -y[,2] * (logfw(eta) - log(sigma)) - (1 - y[,2]) * logSw(eta)
    }

    risk <- function(y, f, w = 1)
        sum(w * loss(y = y, f = f, sigma = sigma))

    ngradient <- function(y, f, w = 1) {
        eta <- (log(y[,1]) - f)/sigma
        nom <- (exp(-eta) + 1)
        ngr <- (y[,2] * (2/nom - 1) + (1 - y[,2])/nom)/sigma
        ngr <- stabilize_ngradient(ngr, w = w, stabilization)
        return(ngr)
    }
    offset <- function(y, w){
        if (!is.null(mu)){
            RET <- mu
        } else {
            if (is.null(sigma))
                sigma <<- 1
            RET <- optimize(risk, interval = c(log(1e-10), log(1e10)),
                            y = y, w = w)$minimum
        }
        return(RET)
    }


    Family(ngradient = ngradient, risk = risk, offset = offset, loss = loss,
           response = function(f) f,
           check_y = function(y) {
               if (!inherits(y, "Surv"))
                   stop("response is not an object of class ", sQuote("Surv"),
                        " but ", sQuote("families = LogLogLSS()"))
               y
           }, name = "Log-Logistic AFT Model: mu (id link)")
}

LogLogSigma <- function (mu = NULL, sigma = NULL, stabilization){
    loss <- function(mu, y, f) {
        logfw <- function(pred)
            dlogis(pred, log = TRUE)
            #exp(pred)/(1 + exp(pred))^2
        logSw <- function(pred)
            pnorm(pred, lower.tail = FALSE, log.p = TRUE)
            #1/(1 + exp(pred))
        eta <- (log(y[,1]) - mu)/exp(f)
        -y[,2] * (logfw(eta) - f) - (1 - y[,2]) * logSw(eta)
    }
    risk <- function(y, f, w = 1)
        sum(w * loss(y = y, f = f, mu = mu))
    ngradient <- function(y, f, w = 1) {
        eta <- (log(y[,1]) - mu)/exp(f)
        ngr <- -(y[,2] + y[,2]*eta -(y[,2]+1)*eta/(1+exp(-eta)))
        ngr <- stabilize_ngradient(ngr, w = w, stabilization)
        return(ngr)
    }
    offset <- function(y, w){
        if (!is.null(sigma)){
            RET <- log(sigma)
        } else {
            if (is.null(mu))
                mu <<- 0
            ## look for starting value of f = log(sigma) in "interval"
            ## i.e. sigma possibly ranges from 1e-10 to 1e10
            RET <- optimize(risk, interval = c(log(1e-10), log(1e10)), y = y, w = w)$minimum
        }
        return(RET)
    }

    Family(ngradient = ngradient, risk = risk, offset = offset, loss = loss,
           response = function(f) exp(f),
           check_y = function(y) {
               if (!inherits(y, "Surv"))
                   stop("response is not an object of class ", sQuote("Surv"),
                        " but ", sQuote("families = LogLogLSS()"))
               y
           },  name = "Log-Logistic AFT Model: sigma (log link)")
}

LogLogLSS <- function(mu = NULL, sigma = NULL,
                      stabilization = c("none", "MAD")) {
    if ((!is.null(sigma) && sigma <= 0))
        stop(sQuote("sigma"), " must be greater than zero")
    stabilization <- check_stabilization(stabilization)
    Families(mu = LogLogMu(mu = mu, sigma = sigma, stabilization = stabilization),
             sigma = LogLogSigma(mu = mu, sigma = sigma, stabilization = stabilization),
#             qfun = qLogLog,
             name = "Log-Log")
}

# nach: qllog (package FAdist)
#qLogLog <- function(p, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
#    exp(qlogis(p = p, location = mu, scale = sigma,
#               lower.tail = lower.tail, log.p = log.p))
#}
# andere Ergebnisse als qllogis (package actuar)


###
# Weibull LSS Family
WeibullMu <- function (mu = NULL, sigma = NULL, stabilization){
    loss <- function(sigma, y, f) {
        logfw <- function(pred)
            pred - exp(pred)
        logSw <- function(pred)
            -exp(pred)
        eta <- (log(y[,1]) - f)/sigma
        -y[,2] * (logfw(eta) -log(sigma)) - (1 - y[,2]) * logSw(eta)
    }
    risk <- function(y, f, w = 1)
        sum(w * loss(y = y, f = f, sigma = sigma))
    ngradient <- function(y, f, w = 1){
        eta <- (log(y[,1]) - f)/sigma
        ngr <- (y[,2] * (exp(eta) - 1) + (1 - y[,2]) * exp(eta))/sigma
        ngr <- stabilize_ngradient(ngr, w = w, stabilization)
        return(ngr)
    }
    offset <- function(y, w){
        if (!is.null(mu)){
            RET <- log(mu)
        } else {
            if (is.null(sigma))
                sigma <<- 1
            RET <- optimize(risk, interval = c(0, max(y[,1], na.rm = TRUE)),
                            y = y, w = w)$minimum
        }
        return(RET)
    }

    Family(ngradient = ngradient, risk = risk, offset = offset, loss = loss,
           response = function(f) f,
           check_y = function(y) {
               if (!inherits(y, "Surv"))
                   stop("response is not an object of class ", sQuote("Surv"),
                        " but ", sQuote("families = WeibullLSS()"))
               y
           }, name = "Weibull AFT Model: mu (id link)")
}


WeibullSigma <- function (mu = NULL, sigma = NULL, stabilization){
    loss <- function(mu, y, f) {
        logfw <- function(pred)
            pred - exp(pred)
        logSw <- function(pred)
            -exp(pred)
        eta <- (log(y[,1]) - mu)/exp(f)
        -y[,2] * (logfw(eta) - f) - (1 - y[,2]) * logSw(eta)
        }
    risk <- function(y, f, w = 1)
        sum(w * loss(y = y, f = f, mu = mu))
    ngradient <- function(y, f, w = 1) {
        eta <- (log(y[,1]) - mu)/exp(f)
        ngr <- -(y[,2] * (eta + 1) - eta * exp(eta))
        ngr <- stabilize_ngradient(ngr, w = w, stabilization)
        return(ngr)
    }
    offset <- function(y, w){
        if (!is.null(sigma)){
            RET <- log(sigma)
        } else {
            if (is.null(mu))
                mu <<- 0
            ## look for starting value of f = log(sigma) in "interval"
            ## i.e. sigma possibly ranges from 1e-10 to 1e10
            RET <- optimize(risk, interval = c(log(1e-10), log(1e10)),
                            y = y, w = w)$minimum
        }
        return(RET)
    }

    Family(ngradient = ngradient, risk = risk, offset = offset, loss = loss,
           response = function(f) exp(f),
           check_y = function(y) {
               if (!inherits(y, "Surv"))
                   stop("response is not an object of class ", sQuote("Surv"),
                        " but ", sQuote("families = WeibullLSS()"))
               y
           }, name = "Weibull AFT Model: sigma (log link)")
}

WeibullLSS <- function(mu = NULL, sigma = NULL,
                       stabilization = c("none", "MAD")) {
    if ((!is.null(sigma) && sigma <= 0))
        stop(sQuote("sigma"), " must be greater than zero")
    stabilization <- check_stabilization(stabilization)
    Families(mu = WeibullMu(mu = mu, sigma = sigma, stabilization = stabilization),
             sigma = WeibullSigma(mu = mu, sigma = sigma, stabilization = stabilization),
             qfun = qWeibull,
             name = "Weibull")
}

qWeibull <- function(p, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
    qweibull(p, scale = mu, shape = sigma,
             lower.tail = lower.tail, log.p = log.p)
}


GaussianMu  <- function(mu = NULL, sigma = NULL, stabilization){

    loss <- function(sigma, y, f) -dnorm(x=y, mean=f, sd=sigma, log=TRUE)
    risk <- function(y, f, w = 1) {
        sum(w * loss(y = y, f = f, sigma = sigma))
    }
    ngradient <- function(y, f, w = 1) {
        ngr <- (1/sigma^2)*(y - f)
        ngr <- stabilize_ngradient(ngr, w = w, stabilization)
        return(ngr)
    }

    offset <- function(y, w){
        if (!is.null(mu)){
            RET <- mu
        } else {
            RET <- weighted.mean(y, w = w, na.rm=TRUE)
        }
        return(RET)
    }

    Family(ngradient = ngradient, risk = risk, loss = loss,
           response = function(f) f, offset=offset,
           name = "Normal distribution: mu(id link)")
}

GaussianSigma  <- function(mu = NULL, sigma = NULL, stabilization){

    loss <-  function(y, f, mu) - dnorm(x=y, mean=mu, sd=exp(f), log=TRUE)
    risk <- function(y, f, w = 1) {
        sum(w * loss(y = y, f = f, mu = mu))
    }
    ngradient <- function(y, f, w = 1) {
        ngr <- (- 1 + exp(-2*f)*((y - mu)^2))
        ngr <- stabilize_ngradient(ngr, w = w, stabilization)
        return(ngr)
    }
    offset <- function(y, w){
        if (!is.null(sigma)){
            RET <- log(sigma)
        } else {
            RET <-  log(weighted.sd(y, w = w, na.rm=TRUE))
        }
        return(RET)
    }

    Family(ngradient = ngradient, risk = risk, loss = loss,
           response = function(f) exp(f), offset=offset,
           name = "Normal distribution: sigma (log link)")
}


GaussianLSS <- function(mu = NULL, sigma = NULL,
                        stabilization = c("none", "MAD")) {
    if ((!is.null(sigma) && sigma <= 0))
        stop(sQuote("sigma"), " must be greater than zero")
    stabilization <- check_stabilization(stabilization)
    Families(mu = GaussianMu(mu=mu, sigma=sigma, stabilization = stabilization),
             sigma = GaussianSigma(mu=mu, sigma=sigma, stabilization = stabilization),
             qfun = qNormal,
             name = "Gaussian")
}


qNormal <- function(p, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
    qnorm(p = p, mean = mu, sd = sigma, lower.tail = lower.tail, log.p = log.p)
}

GammaMu <-function (mu = NULL, sigma = NULL, stabilization) {
    loss <-  function(sigma, y, f) {
        lgamma(sigma) + sigma * y * exp(-f) - sigma * log(y) -
            sigma * log(sigma) + sigma * f + log(y)
    }
    risk <- function(y, f, w = 1) {
        sum(w * loss(y = y, f = f, sigma = sigma))
    }
    ngradient <- function(y, f, w = 1) {
        ngr <- sigma * y * exp(-f) - sigma
        ngr <- stabilize_ngradient(ngr, w = w, stabilization)
        return(ngr)
    }
    offset <- function(y, w) {
        if (!is.null(mu)) {
            RET <- log(mu)
        }
        else {
            if (is.null(sigma))
                sigma <<- mean(y)^2/(var(y))
            RET <- optimize(risk, interval = c(0, max(y^2, na.rm = TRUE)),
            y = y, w = w)$minimum
        }
        return(RET)
    }
    check_y <- function(y) {
        if (!is.numeric(y) || !is.null(dim(y)))
            stop("response is not a numeric vector but ", sQuote("GammaLSS()"))
        if (any(y < 0))
            stop("response is not positive but ", sQuote("GammaLSS()"))
        y
    }
    Family(ngradient = ngradient, risk = risk, loss = loss,
           response = function(f) exp(f), offset = offset, check_y = check_y,
           name = "Gamma distribution: mu(log link)")
}

GammaSigma <- function(mu = NULL, sigma = NULL, stabilization) {
    loss <-  function(mu, y, f) {
        lgamma(exp(f)) + (exp(f) * y )/mu - exp(f) * log(y) -
            f * exp(f) + exp(f) * log(mu) + log(y)
    }
    risk <- function(y, f, w = 1){
        sum(w * loss(y = y, f = f, mu = mu))
    }
    ngradient <- function(y, f, w = 1) {
        ngr <- - digamma(exp(f))*exp(f) + (f+1)*exp(f) - log(mu)*exp(f) +
            exp(f)*log(y) - (y*exp(f))/mu
        ngr <- stabilize_ngradient(ngr, w = w, stabilization)
        return(ngr)
    }
    offset <- function(y, w) {
        if (!is.null(sigma)) {
            RET <- log(sigma)
        }
        else {
            if (is.null(mu))
                mu <<- mean(y)
            RET <- optimize(risk, interval = c(0, max(y, na.rm = TRUE)),
            y = y, w = w)$minimum
        }
        return(RET)
    }
    check_y <- function(y) {
        if (!is.numeric(y) || !is.null(dim(y)))
            stop("response is not a numeric vector but ", sQuote("GammaLSS()"))
        if (any(y < 0))
            stop("response is not positive but ", sQuote("GammaLSS()"))
        y
    }
    Family(ngradient = ngradient, risk = risk, loss = loss,
           response = function(f) exp(f), offset = offset, check_y = check_y,
           name = "Gamma distribution: sigma(log link)")
}

GammaLSS <- function (mu = NULL, sigma = NULL,
                      stabilization = c("none", "MAD")) {
    if ((!is.null(sigma) && sigma <= 0))
        stop(sQuote("sigma"), " must be greater than zero")
    if ((!is.null(mu) && mu <= 0))
        stop(sQuote("mu"), " must be greater than zero")
    stabilization <- check_stabilization(stabilization)
    Families(mu = GammaMu(mu = mu, sigma = sigma, stabilization = stabilization),
             sigma = GammaSigma(mu = mu, sigma = sigma, stabilization = stabilization),
             qfun = qGamma,
             name = "Gamma")
}

## we use almost the same parameterization as GA but with sigma = sqrt(1/sigma)
qGamma <- function(p, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
    ## require gamlss.dist
    if (!requireNamespace("gamlss.dist", quietly = TRUE))
        stop("Please install package 'gamlss.dist' for using qGamma.")
    gamlss.dist::qGA(p = p, mu = mu, sigma = sqrt(1/sigma), lower.tail = lower.tail, log.p = log.p)
}



BetaMu <- function(mu = NULL, phi = NULL, stabilization){

    # loss is negative log-Likelihood, f is the parameter to be fitted with
    # logit link -> exp(f) = mu
    loss <- function(phi, y, f) {
        - 1 * (lgamma(phi) - lgamma(plogis(f) * phi) -
              lgamma((1 - plogis(f)) * phi) + (plogis(f) * phi - 1) * log(y) +
              ((1 - plogis(f)) * phi - 1) * log(1 - y))
    }
    # risk is sum of loss
    risk <- function(y, f, w = 1) {
        sum(w * loss(y = y, f = f, phi = phi))
        #sqrt(mean(w * (y - plogis(f))^2))
    }
    # ngradient is the negative derivate w.r.t. mu
    ngradient <- function(y, f, w = 1) {
        ngr <- +1 * exp(f)/(1 + exp(f))^2 * (phi * (qlogis(y) - (digamma(plogis(f) * phi) -
              digamma((1 - plogis(f)) * phi)))) # Nachdifferenzieren? -> nein, da nach mu ableiten und nicht nach beta
        ngr <- stabilize_ngradient(ngr, w = w, stabilization)
        return(ngr)
    }
    offset <- function(y, w) {
        if (!is.null(mu)) {
            RET <- qlogis(mu)
        }
        else {
            if (is.null(phi))
               phi <<-  mean(y) * (1 - mean(y)) /var(y) - 1
            ### look for starting value of f = qlogis(mu) in "interval"
            ### i.e. mu ranges from 0.000001 to 0.999999
            RET <- optimize(risk, interval = c(qlogis(0.000001), qlogis(0.999999)), y = y,
                w = w)$minimum
        }
        return(RET)
    }
    # use the Family constructor of mboost
    Family(ngradient = ngradient, risk = risk, loss = loss, offset = offset,
           response = function(f) plogis(f),
           check_y <- function(y) {
               if (!is.numeric(y) || !is.null(dim(y)))
                   stop("response must be a numeric vector")
               if (any(y <= 0) | any(y >= 1))
                   stop("response must be >0 and <1")
               y
           },
           name = "Beta-distribution: mu (logit link)")
}

# Sub-Family for phi
BetaPhi <- function(mu = NULL, phi = NULL, stabilization){

    # loss is negative log-Likelihood, f is the parameter to be fitted with
    # log-link: exp(f) = phi
    loss <- function(mu, y, f) {
        #y <- (y * (length(y) - 1) + 0.5)/length(y)
        -1 * (lgamma(exp(f)) - lgamma(mu * exp(f)) -
              lgamma((1 - mu) * exp(f)) + (mu * exp(f) - 1) * log(y) +
              ((1 - mu) * exp(f) - 1) * log(1 - y))
    }
    # risk is sum of loss
    risk <- function(y, f, w = 1) {
        sum(w * loss(y = y, f = f, mu = mu))
    }
    # ngradient is the negative derivate w.r.t. phi
    ngradient <- function(y, f, w = 1) {
        #y <- (y * (length(y) - 1) + 0.5)/length(y)
        ngr <- +1 * exp(f) * (mu * (qlogis(y) - (digamma(mu * exp(f)) - digamma((1 - mu) * exp(f)))) +
              digamma(exp(f)) - digamma((1 - mu) * exp(f)) + log(1 - y))
        ngr <- stabilize_ngradient(ngr, w = w, stabilization)
        return(ngr)
    }
    offset <- function(y, w) {
        if (!is.null(phi)) {
            RET <- log(phi)
        }
        else {
            if (is.null(mu))
                mu <<- mean(y)
            ### look for starting value of f = log(phi) in "interval"
            ### i.e. phi possibly ranges from 1e-10 to 1e10
            RET <- optimize(risk, interval = c(log(1e-10), log(1e10)),
                y = y, w = w)$minimum
        }
        return(RET)
    }
    # use the Family constructor of mboost
    Family(ngradient = ngradient, risk = risk, loss = loss, offset = offset,
           response = function(f) exp(f),
           check_y <- function(y) {
               if (!is.numeric(y) || !is.null(dim(y)))
                   stop("response must be a numeric vector")
               if (any(y <= 0) | any(y >= 1))
                   stop("response must be >0 and <1")
               y
           },
           name = "Beta-distribution: phi (log link)")
}

# families object for new distribution
BetaLSS <- function (mu = NULL, phi = NULL,
                     stabilization = c("none", "MAD")) {
    stabilization <- check_stabilization(stabilization)
    Families(mu = BetaMu(mu = mu, phi = phi, stabilization = stabilization),
             phi = BetaPhi(mu = mu, phi = phi, stabilization = stabilization),
             qfun = qBeta,
             name = "Beta")
}

## we use almost the same parameterization as BE but with sigma = 1/sqrt(phi + 1)
qBeta <- function(p, mu = 0, phi = 1, lower.tail = TRUE, log.p = FALSE) {
    ## require gamlss.dist
    if (!requireNamespace("gamlss.dist", quietly = TRUE))
        stop("Please install package 'gamlss.dist' for using qBeta.")
    gamlss.dist::qBE(p = p, mu = mu, sigma = 1/sqrt(phi + 1), lower.tail = lower.tail, log.p = log.p)
}

# Zero-inflated Poisson model
ZIPoLSS <- function(mu = NULL, sigma = NULL,
                    stabilization = c("none", "MAD")) {
    fam <- as.families(fname = "ZIP", mu = mu, sigma = sigma, stabilization = stabilization)

    fam$mu@name <- "Zero-inflated Poisson model: count data component"
    fam$sigma@name <- "Zero-inflated Poisson model: zero component"

    fam
}

# Zero-inflated negative binomial model
ZINBLSS <- function(mu = NULL, sigma = NULL, nu = NULL,
                    stabilization = c("none", "MAD")) {
    fam <- as.families(fname = "ZINBI", mu = mu, sigma = sigma, nu = nu, stabilization = stabilization)

    fam$mu@name <- "Zero-inflated negative binomial model: location parameter for count data component"
    fam$sigma@name <- "Zero-inflated negative binomial model: scale parameter for count data component"
    fam$nu@name <- "Zero-inflated negative binomial model: zero component"

    fam
}
