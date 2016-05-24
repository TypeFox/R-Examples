dskewhyp <- function (x, mu = 0, delta = 1, beta = 1, nu = 1,
                      param = c(mu,delta,beta,nu), log = FALSE,
                      tolerance = .Machine$double.eps^0.5)
{
    parResult <- skewhypCheckPars(param)
    case <- parResult$case
    errMessage <- parResult$errMessage
    if (case == "error")
        stop(errMessage)
    param <- as.numeric(param)
    mu <- param[1]
    delta <- param[2]
    beta <- param[3]
    nu <- param[4]
    if (abs(beta) > tolerance) {
        ldskewhyp <- ((1 - nu)/2) * log(2) + nu * log(delta) +
                     ((nu + 1)/2) * log(abs(beta)) +
                     log(besselK(x = sqrt(beta^2*(delta^2 + (x - mu)^2)),
                                 nu = (nu + 1)/2, expon.scaled = TRUE)) -
                     sqrt(beta^2 * (delta^2 + (x - mu)^2)) +
                     beta * (x - mu) - lgamma(nu/2) - log(pi)/2 -
                     ((nu + 1)/2) * log(delta^2 + (x - mu)^2)/2
    }
    else {
        ldskewhyp <- lgamma((nu + 1)/2) - log(pi)/2 - log(delta) -
            lgamma(nu/2) - ((nu + 1)/2) * log(1 + ((x - mu)^2)/delta^2)
    }
    if (log == TRUE) {
        return(ldskewhyp)
    }
    else {
        return(exp(ldskewhyp))
    }
}


pskewhyp <- function (q, mu = 0, delta = 1, beta = 1, nu = 1,
                      param = c(mu, delta, beta, nu), log.p = FALSE,
                      lower.tail = TRUE, subdivisions = 100,
                      intTol = .Machine$double.eps^0.25,
                      valueOnly = TRUE, ...)
{
    parResult <- skewhypCheckPars(param)
    case <- parResult$case
    errMessage <- parResult$errMessage
    if (case == "error")
        stop(errMessage)
    mu <- param[1]
    delta <- param[2]
    beta <- param[3]
    nu <- param[4]


    distMode <- skewhypMode(param = param)
    qLess <- which((q <= distMode)&(is.finite(q)))
    qGreater <- which((q > distMode)&(is.finite(q)))
    prob <- rep(NA, length(q))
    err <- rep(NA, length(q))

    prob[q == -Inf] <- 0
    ## looks wrong but will be changed later: upper tail prob
    prob[q == Inf] <- 0
    err[q %in% c(-Inf,Inf)] <- 0

    dskewhypInt <- function(q) {
        dskewhyp(q, param = param)
    }

    for (i in qLess) {
        intRes <- integrate(dskewhypInt, -Inf, q[i],
                          subdivisions = subdivisions,
                          rel.tol = intTol, ...)
        prob[i] <- intRes$value
        err[i] <- intRes$abs.error
    }
    for (i in qGreater) {
        intRes <- integrate(dskewhypInt, q[i], Inf,
                            subdivisions = subdivisions,
                            rel.tol = intTol, ...)
        prob[i] <- intRes$value
        err[i] <- intRes$abs.error
    }

    if (lower.tail == TRUE) {
        if (length(q > distMode) > 0){
            ##cat("q =", q, "distMode = ", distMode, "prob = ", prob, "\n")
            prob[q > distMode] <- 1 - prob[q > distMode]
        }
    } else {
        if (length(q <= distMode) > 0){
            prob[q <= distMode] <- 1 - prob[q <= distMode]
        }
    }

    if (log.p == TRUE) {
        prob <- log(prob)
    }

    ifelse(valueOnly, return(prob),
           return(list(value = prob, error = err)))
}

qskewhyp <- function (p, mu = 0, delta = 1, beta = 1, nu = 1,
                      param = c(mu,delta, beta, nu),
                      lower.tail = TRUE, log.p = FALSE,
                      method = c("spline","integrate"),
                      nInterpol = 501, uniTol = .Machine$double.eps^0.25,
                      subdivisions = 100, intTol = uniTol, ...)
{
    parResult <- skewhypCheckPars(param)
    case <- parResult$case
    errMessage <- parResult$errMessage
    if (case == "error")
        stop(errMessage)
    if(!lower.tail){
      p <- 1 - p
      lower.tail == TRUE
    }
    method <- match.arg(method)
    param <- as.numeric(param)
    mu <- param[1]
    delta <- param[2]
    beta <- param[3]
    nu <- param[4]
    distMode <- skewhypMode(param = param)
    pModeDist <- pskewhyp(distMode, param = param,
                          intTol = intTol)
    xRange <- skewhypCalcRange(param = param, tol = 10^(-7))

    quant <- rep(NA, length(p))
    invalid <- which(p < 0 | p > 1)
    pFinite <- which((p > 0) & (p < 1))

    ## deal with limit values
    quant[p == 0] <- -Inf
    quant[p == 1] <- Inf

    if (method == "integrate"){
        less <-  which((p <= pModeDist) & (p > 0))
        if (length(less) > 0){
            pLow <- min(p[less])
            xLow <- distMode  -
                skewhypStepSize(delta, delta, beta, nu, "left")
            while (pskewhyp(xLow, param = param) >= pLow){
                xLow <- xLow -
                    skewhypStepSize(distMode - xLow, delta, beta, nu, "left")
            }
            xRange <- c(xLow,distMode)
            zeroFn <- function(x, param, p) {
                return(pskewhyp(x, param = param,
                                subdivisions = subdivisions,
                                intTol = intTol) - p)
            }
            for (i in less){
                quant[i] <- uniroot(zeroFn, param = param, p = p[i],
                                    interval = xRange, tol = uniTol)$root
            }
        }
        greater <-  which((p > pModeDist) & (p < 1))
        p[greater] <- 1 - p[greater]
        if (length(greater) > 0){
            pHigh <- min(p[greater])
            xHigh <- distMode  +
                skewhypStepSize(delta, delta, beta, nu, "right")
            while(pskewhyp(xHigh, param = param, lower.tail = FALSE)
                  >= pHigh){
                xHigh <- xHigh +
                    skewhypStepSize(xHigh - distMode, delta, beta, nu, "right")
            }
            xRange <- c(distMode,xHigh)
            zeroFn <- function(x, param, p) {
                return(pskewhyp(x, param = param, lower.tail = FALSE,
                                subdivisions = subdivisions,
                                intTol = intTol) - p)
            }
            for (i in greater){
                quant[i] <- uniroot(zeroFn, param = param, p = p[i],
                                    interval = xRange, tol = uniTol)$root
            }
        }
    } else if (method == "spline"){
        inRange <- which((p > pskewhyp(xRange[1], param = param)) &
                         (p < pskewhyp(xRange[2], param = param)))
        small <- which((p <= pskewhyp(xRange[1], param = param)) &
                       (p > 0))
        large <- which((p >= pskewhyp(xRange[2], param = param)) &
                       (p < 1))
        extreme <- c(small,large)
        xVals <- seq(xRange[1], xRange[2], length.out = nInterpol)
        yVals <- pskewhyp(xVals, param = param, subdivisions = subdivisions,
                          intTol = intTol)
        splineFit <- splinefun(xVals, yVals)
        zeroFn <- function(x, p) {
            return(splineFit(x) - p)
        }
        for (i in inRange){
            quant[i] <- uniroot(zeroFn, p = p[i],
                                interval = xRange, tol = uniTol)$root
        }

        if (length(extreme) > 0){
            quant[extreme] <- qskewhyp(p[extreme], param = param,
                                       lower.tail = lower.tail, log.p = log.p,
                                       method = "integrate",
                                       nInterpol = nInterpol, uniTol = uniTol,
                                       subdivisions = subdivisions,
                                       intTol = intTol, ...)
        }

    }

    return(quant)
}

rskewhyp <- function (n, mu = 0, delta = 1, beta = 1, nu = 1,
                      param = c(mu,delta,beta,nu), log = FALSE)
{
    parResult <- skewhypCheckPars(param)
    case <- parResult$case
    errMessage <- parResult$errMessage
    if (case == "error")
        stop(errMessage)
    param <- as.numeric(param)
    mu <- param[1]
    delta <- param[2]
    beta <- param[3]
    nu <- param[4]
    if (log == TRUE)
        stop("This function is not yet implemented")
    y <- 1/rgamma(n, shape = nu/2, scale = 2/delta^2)
    sigma <- sqrt(y)
    z <- rnorm(n)
    rskewhyp <- mu + beta * sigma^2 + sigma * z
    return(rskewhyp)
}

ddskewhyp <- function (x, mu = 0, delta = 1, beta = 1, nu = 1,
                       param = c(mu,delta,beta,nu), log = FALSE,
                       tolerance = .Machine$double.eps^0.5)
{
    parResult <- skewhypCheckPars(param)
    case <- parResult$case
    errMessage <- parResult$errMessage
    if (case == "error")
        stop(errMessage)
    param <- as.numeric(param)
    mu <- param[1]
    delta <- param[2]
    beta <- param[3]
    nu <- param[4]
    if (log == TRUE)
        stop("This function is not yet implemented")
    if (abs(beta) > tolerance) {
        ddskewhyp <- 1/2 * 2^(1/2 - 1/2 * nu) * delta^nu * abs(beta)^(1/2 *
            nu + 1/2) * (-besselK(nu = 1/2 * nu + 3/2, x = (beta^2 *
            (delta^2 + (x - mu)^2))^(1/2)) + (1/2 * nu + 1/2)/(beta^2 *
            (delta^2 + (x - mu)^2))^(1/2) * besselK(nu = 1/2 *
            nu + 1/2, x = (beta^2 * (delta^2 + (x - mu)^2))^(1/2)))/(beta^2 *
            (delta^2 + (x - mu)^2))^(1/2) * beta^2 * (2 * x -
            2 * mu) * exp(beta * (x - mu))/gamma(1/2 * nu)/pi^(1/2)/(((delta^2 +
            x^2 - 2 * x * mu + mu^2)^(1/2))^(1/2 * nu + 1/2)) +
            2^(1/2 - 1/2 * nu) * delta^nu * abs(beta)^(1/2 *
                nu + 1/2) * besselK(nu = 1/2 * nu + 1/2, x = (beta^2 *
                (delta^2 + (x - mu)^2))^(1/2)) * beta * exp(beta *
                (x - mu))/gamma(1/2 * nu)/pi^(1/2)/(((delta^2 +
                x^2 - 2 * x * mu + mu^2)^(1/2))^(1/2 * nu + 1/2)) -
            1/2 * 2^(1/2 - 1/2 * nu) * delta^nu * abs(beta)^(1/2 *
                nu + 1/2) * besselK(nu = 1/2 * nu + 1/2, x = (beta^2 *
                (delta^2 + (x - mu)^2))^(1/2)) * exp(beta * (x -
                mu))/gamma(1/2 * nu)/pi^(1/2)/(((delta^2 + x^2 -
                2 * x * mu + mu^2)^(1/2))^(1/2 * nu + 1/2)) *
                (1/2 * nu + 1/2)/(delta^2 + x^2 - 2 * x * mu +
                mu^2) * (2 * x - 2 * mu)
    }
    else {
        ddskewhyp <- 2 * gamma(1/2 * nu + 1/2)/pi^(1/2)/delta^3/gamma(1/2 *
            nu) * (1 + (x - mu)^2/delta^2)^(-1/2 * nu - 1/2) *
            (-1/2 * nu - 1/2) * (x - mu)/(1 + (x - mu)^2/delta^2)
    }
    return(ddskewhyp)
}
