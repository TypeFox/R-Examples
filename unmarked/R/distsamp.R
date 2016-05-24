
distsamp <- function(formula, data,
    keyfun=c("halfnorm", "exp", "hazard", "uniform"),
    output=c("density", "abund"), unitsOut=c("ha", "kmsq"), starts=NULL,
    method="BFGS", se = TRUE, engine = c("C", "R"),
    rel.tol=0.001, ...)
{
    engine <- match.arg(engine)
    keyfun <- match.arg(keyfun)
#    if(engine=="C" && !(keyfun %in% c("halfnorm", "exp", "uniform"))) {
#        engine <- "R"
#        warning("C engine not available for hazard model, using R instead")
#    }
    output <- match.arg(output)
    unitsOut <- match.arg(unitsOut)
    db <- data@dist.breaks
    tlength <- data@tlength
    survey <- data@survey
    w <- diff(db)
    unitsIn <- data@unitsIn
    designMats <- getDesign(data, formula)
    X <- designMats$X; V <- designMats$V; y <- designMats$y
    X.offset <- designMats$X.offset; V.offset <- designMats$V.offset
    if(is.null(X.offset))
        X.offset <- rep(0, nrow(X))
    if(is.null(V.offset))
        V.offset <- rep(0, nrow(V))
    M <- nrow(y)
    J <- ncol(y)

    u <- a <- matrix(NA, M, J)
    switch(survey,
        line = {
            for(i in 1:M) {
                a[i,] <- tlength[i] * w
                u[i,] <- a[i,] / sum(a[i,])
                }
            },
        point = {
            for(i in 1:M) {
                a[i, 1] <- pi*db[2]^2
                for(j in 2:J)
                    a[i, j] <- pi*db[j+1]^2 - sum(a[i, 1:(j-1)])
                u[i,] <- a[i,] / sum(a[i,])
                }
            })
    switch(survey,
        line = A <- rowSums(a) * 2,
        point = A <- rowSums(a))
    switch(unitsIn,
        m = A <- A / 1e6,
        km = A <- A)
    switch(unitsOut,
        ha = A <- A * 100,
        kmsq = A <- A)

    lamParms <- colnames(X)
    detParms <- colnames(V)
    nAP <- length(lamParms)
    nDP <- length(detParms)
    nP <- nAP + nDP
    cp <- matrix(NA, M, J)
    switch(keyfun,
           halfnorm = {
               altdetParms <- paste("sigma", colnames(V), sep="")
               if(is.null(starts)) {
                   starts <- c(rep(0, nAP), log(max(db)), rep(0, nDP-1))
                   names(starts) <- c(lamParms, detParms)
               } else {
                   if(is.null(names(starts))) names(starts) <-
                       c(lamParms, detParms)
               }
           },
           exp = {
               altdetParms <- paste("rate", colnames(V), sep="")
               if(is.null(starts)) {
                   starts <- c(rep(0, nAP), 0, rep(0, nDP-1))
                   names(starts) <- c(lamParms, detParms)
               } else {
                   if(is.null(names(starts))) names(starts) <-
                       c(lamParms, detParms)
               }
           },
           hazard = {
               nDP <- length(detParms)
               nP <- nAP + nDP + 1
               altdetParms <- paste("shape", colnames(V), sep="")
               if(is.null(starts)) {
                   starts <- c(rep(0, nAP), log(median(db)),
                               rep(0, nDP-1), 1)
                   names(starts) <- c(lamParms, detParms, "scale")
               } else {
                   if(is.null(names(starts)))
                       names(starts) <- c(lamParms, detParms, "scale")
               }
           },
           uniform = {
               detParms <- character(0)
               altdetParms <- character(0)
               nDP <- 0
               if(is.null(starts)) {
                   starts <- rep(0, length(lamParms))
                   names(starts) <- lamParms
               } else {
                   if(is.null(names(starts))) names(starts) <- lamParms
               }
           })

    if(engine=="R") {
    switch(keyfun,
    halfnorm = {
        nll <- function(param) {
            sigma <- drop(exp(V %*% param[(nAP+1):nP] + V.offset))
            lambda <- drop(exp(X %*% param[1:nAP] + X.offset))
            if(identical(output, "density"))
                lambda <- lambda * A
            for(i in 1:M) {
                switch(survey,
                line = {
                    f.0 <- 2 * dnorm(0, 0, sd=sigma[i])
                    int <- 2 * (pnorm(db[-1], 0, sd=sigma[i]) -
                        pnorm(db[-(J+1)], 0, sd=sigma[i]))
                    cp[i,] <- int / f.0 / w
                    },
                point = {
                    for(j in 1:J) {
                        int <- integrate(grhn, db[j], db[j+1], sigma=sigma[i],
                            stop.on.error=FALSE)
                        mess <- int$message
                        if(identical(mess, "OK"))
                            cp[i, j] <- int$value * 2*pi / a[i,j]
                        else {
                            cp[i, j] <- NA
                            }
                        }
                    })
                cp[i,] <- cp[i,] * u[i,]
                }
            ll <- dpois(y, lambda * cp, log=TRUE)
            -sum(ll)
            }},
    exp = {
        nll <- function(param) {
            rate <- drop(exp(V %*% param[(nAP+1):nP] + V.offset))
            lambda <- drop(exp(X %*% param[1:nAP] + X.offset))
            if(identical(output, "density"))
                lambda <- lambda * A
            for(i in 1:M) {
                switch(survey,
                line = {
                    for(j in 1:J) {
                        int <- integrate(gxexp, db[j], db[j+1], rate=rate[i],
                            stop.on.error=FALSE)
                        mess <- int$message
                        if(identical(mess, "OK"))
                            cp[i, j] <- int$value / w[j]
                        else {
                            cp[i, j] <- NA
                            }
                        }},
                point = {
                    for(j in 1:J) {
                        int <- integrate(grexp, db[j], db[j+1], rate=rate[i],
                            stop.on.error=FALSE)
                        mess <- int$message
                        if(identical(mess, "OK"))
                            cp[i, j] <- int$value * 2*pi / a[i,j]
                        else {
                            cp[i, j] <- NA
                            }
                        }
                    })
                cp[i,] <- cp[i,] * u[i,]
                }
            ll <- dpois(y, lambda * cp, log=TRUE)
            -sum(ll)
            }},
    hazard = {
        nll <- function(param) {
            shape <- drop(exp(V %*% param[(nAP+1):(nP-1)] + V.offset))
            scale <- drop(exp(param[nP]))
            lambda <- drop(exp(X %*% param[1:nAP] + X.offset))
            if(identical(output, "density"))
                lambda <- lambda * A
            for(i in 1:M) {
                switch(survey,
                line = {
                    for(j in 1:J) {
                        int <- integrate(gxhaz, db[j], db[j+1], shape=shape[i],
                            scale=scale, stop.on.error=FALSE)
                        mess <- int$message
                        if(identical(mess, "OK"))
                            cp[i, j] <- int$value / w[j]
                        else {
                            cp[i, j] <- NA
                            }
                        }},
                point = {
                    for(j in 1:J) {
                        int <- integrate(grhaz, db[j], db[j+1], shape=shape[i],
                            scale=scale, stop.on.error=FALSE)
                        mess <- int$message
                        if(identical(mess, "OK"))
                            cp[i, j] <- int$value * 2*pi / a[i,j]
                        else {
                            cp[i, j] <- NA
                            }
                        }

                    })
                cp[i,] <- cp[i,] * u[i,]
                }
            ll <- dpois(y, lambda * cp, log=TRUE)
            -sum(ll)
            }},
    uniform = {
        nll <- function(param) {
            lambda <- drop(exp(X %*% param + X.offset))
            if(identical(output, "density"))
                lambda <- lambda * A
            ll <- dpois(y, lambda * u, log=TRUE)
            -sum(ll)
            }
        })
    } else if(engine=="C") {
        nll <- function(param) {
            beta.lam <- param[1:nAP]
            if(identical(keyfun, "hazard")) {
                beta.sig <- param[(nAP+1):(nP-1)]
                scale <- exp(param[nP])
            } else {
                beta.sig <- param[(nAP+1):nP]
                scale <- -99.0
            }
            lambda <- drop(exp(X %*% beta.lam + X.offset))
            if(identical(output, "density"))
                lambda <- lambda * A
            sigma <- drop(exp(V %*% beta.sig + V.offset))
            .Call("nll_distsamp",
                  y, lambda, sigma, scale,
                  a, u, w, db,
                  keyfun, survey, rel.tol,
                  PACKAGE="unmarked")
        }
    }
    fm <- optim(starts, nll, method=method, hessian=se, ...)
    opt <- fm
    ests <- fm$par
    if(se) {
        covMat <- tryCatch(solve(fm$hessian), error=function(x)
        simpleError("Hessian is singular. Try using fewer covariates or providing starting values."))
        if(class(covMat)[1] == "simpleError") {
            print(covMat$message)
            covMat <- matrix(NA, nP, nP)
            }
        }
    else
        covMat <- matrix(NA, nP, nP)
    estsAP <- ests[1:nAP]
    if(keyfun == "hazard") {
        estsDP <- ests[(nAP+1):(nP-1)]
        estsScale <- ests[nP]
        }
    else
        estsDP <- ests[(nAP+1):nP]
    covMatAP <- covMat[1:nAP, 1:nAP, drop=F]
    if(keyfun=="hazard") {
        covMatDP <- covMat[(nAP+1):(nP-1), (nAP+1):(nP-1), drop=F]
        covMatScale <- covMat[nP, nP, drop=F]
        }
    else if(keyfun!="uniform")
        covMatDP <- covMat[(nAP+1):nP, (nAP+1):nP, drop=F]
    names(estsDP) <- altdetParms
    fmAIC <- 2 * fm$value + 2 * nP
    stateName <- switch(output, abund = "Abundance", density = "Density")
    stateEstimates <- unmarkedEstimate(name = stateName,
        short.name = "lam", estimates = estsAP, covMat = covMatAP,
        invlink = "exp", invlinkGrad = "exp")
    if(keyfun != "uniform") {
        detEstimates <- unmarkedEstimate(name = "Detection", short.name = "p",
        estimates = estsDP, covMat = covMatDP, invlink = "exp",
        invlinkGrad = "exp")
        if(keyfun != "hazard")
            estimateList <- unmarkedEstimateList(list(
                state=stateEstimates, det=detEstimates))
        else {
            scaleEstimates <- unmarkedEstimate(name = "Hazard-rate(scale)",
                short.name = "p", estimates = estsScale,
                covMat = covMatScale, invlink = "exp", invlinkGrad = "exp")
            estimateList <- unmarkedEstimateList(list(state=stateEstimates,
                det=detEstimates, scale=scaleEstimates))
            }
        }
    else
        estimateList <- unmarkedEstimateList(list(state=stateEstimates))
    dsfit <- new("unmarkedFitDS", fitType = "distsamp", call = match.call(),
        opt = opt, formula = formula, data = data, keyfun=keyfun,
        sitesRemoved = designMats$removed.sites, unitsOut=unitsOut,
        estimates = estimateList, AIC = fmAIC, negLogLike = fm$value,
        nllFun = nll, output=output)
    return(dsfit)
}


# Detection functions

gxhn <- function(x, sigma) exp(-x^2/(2 * sigma^2))
gxexp <- function(x, rate) exp(-x / rate)
gxhaz <- function(x, shape, scale)  1 - exp(-(x/shape)^-scale)
grhn <- function(r, sigma) exp(-r^2/(2 * sigma^2)) * r
grexp <- function(r, rate) exp(-r / rate) * r
grhaz <- function(r, shape, scale)  (1 - exp(-(r/shape)^-scale)) * r

dxhn <- function(x, sigma)
	gxhn(x=x, sigma=sigma) / integrate(gxhn, 0, Inf, sigma=sigma)$value
drhn <- function(r, sigma)
	grhn(r=r, sigma=sigma) / integrate(grhn, 0, Inf, sigma=sigma)$value
dxexp <- function(x, rate)
	gxexp(x=x, rate=rate) / integrate(gxexp, 0, Inf, rate=rate)$value
drexp <- function(r, rate)
	grexp(r=r, rate=rate) / integrate(grexp, 0, Inf, rate=rate)$value
dxhaz <- function(x, shape, scale)
	gxhaz(x=x, shape=shape, scale=scale) / integrate(gxhaz, 0, Inf,
		shape=shape, scale=scale)$value
drhaz <- function(r, shape, scale)
	grhaz(r=r, shape=shape, scale=scale) / integrate(grhaz, 0, Inf,
		shape=shape, scale=scale)$value


