
gdistsamp <- function(lambdaformula, phiformula, pformula, data,
    keyfun=c("halfnorm", "exp", "hazard", "uniform"),
    output=c("abund", "density"), unitsOut=c("ha", "kmsq"),
    mixture=c('P', 'NB'), K,
    starts, method = "BFGS", se = TRUE, rel.tol=1e-4,
    ...)
{
if(!is(data, "unmarkedFrameGDS"))
    stop("Data is not of class unmarkedFrameGDS.")

keyfun <- match.arg(keyfun)
if(!keyfun %in% c("halfnorm", "exp", "hazard", "uniform"))
    stop("keyfun must be 'halfnorm', 'exp', 'hazard', or 'uniform'")
output <- match.arg(output)
unitsOut <- match.arg(unitsOut)
db <- data@dist.breaks
w <- diff(db)
tlength <- data@tlength
survey <- data@survey
unitsIn <- data@unitsIn
mixture <- match.arg(mixture)

formlist <- list(lambdaformula = lambdaformula, phiformula = phiformula,
    pformula = pformula)
form <- as.formula(paste(unlist(formlist), collapse=" "))
D <- getDesign(data, formula = form)

Xlam <- D$Xlam
Xphi <- D$Xphi
Xdet <- D$Xdet
y <- D$y  # MxJT

Xlam.offset <- D$Xlam.offset
Xphi.offset <- D$Xphi.offset
Xdet.offset <- D$Xdet.offset
if(is.null(Xlam.offset)) Xlam.offset <- rep(0, nrow(Xlam))
if(is.null(Xphi.offset)) Xphi.offset <- rep(0, nrow(Xphi))
if(is.null(Xdet.offset)) Xdet.offset <- rep(0, nrow(Xdet))

M <- nrow(y)
T <- data@numPrimary
R <- ncol(y)
J <- R / T

y <- array(y, c(M, J, T))
y <- aperm(y, c(1,3,2))
yt <- apply(y, 1:2, function(x) {
    if(all(is.na(x)))
        return(NA)
    else return(sum(x, na.rm=TRUE))
    })

if(missing(K) || is.null(K)) K <- max(yt, na.rm=TRUE) + 100
k <- 0:K
lk <- length(k)

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


lamPars <- colnames(Xlam)
if(T==1) {
    phiPars <- character(0)
    nPP <- 0
    }
else {
    phiPars <- colnames(Xphi)
    nPP <- ncol(Xphi)
    }
if(identical(keyfun, "uniform")) {
    nDP <- 0
    detPars <- character(0)
    }
else {
    nDP <- ncol(Xdet)
    detPars <- colnames(Xdet)
    }
if(identical(keyfun, "hazard")) {
    nSP <- 1
    scalePar <- "scale"
    }
else {
    nSP <- 0
    scalePar <- character(0)
    }
if(identical(mixture, "NB")) {
    nOP <- 1
    nbPar <- "alpha"
    }
else {
    nOP <- 0
    nbPar <- character(0)
    }

nLP <- ncol(Xlam)
nP  <- nLP + nPP + nDP + nSP + nOP

cp <- array(as.numeric(NA), c(M, T, J+1))
g <- matrix(as.numeric(NA), M, lk)

lfac.k <- lgamma(k+1)
kmyt <- array(NA, c(M, T, lk))
lfac.kmyt <- array(0, c(M, T, lk))
fin <- matrix(NA, M, lk)
naflag <- array(NA, c(M, T, J))
for(i in 1:M) {
    fin[i, ] <- k - max(yt[i,], na.rm=TRUE) >= 0
    for(t in 1:T) {
        naflag[i,t,] <- is.na(y[i,t,])
        if(!all(naflag[i,t,])) {
            kmyt[i,t,] <- k - yt[i,t]
            lfac.kmyt[i, t, fin[i,]] <- lgamma(kmyt[i, t, fin[i,]] + 1)
            }
        }
    }

switch(keyfun,
halfnorm = {
    altdetParms <- paste("sigma", colnames(Xdet), sep="")
    if(missing(starts)) {
        starts <- rep(0, nP)
        starts[nLP+nPP+1] <- log(max(db))
        }

    nll <- function(pars) {
        lambda <- exp(Xlam %*% pars[1:nLP] + Xlam.offset)
        if(identical(output, "density"))
            lambda <- lambda * A

        if(T==1)
            phi <- matrix(1, M, T)
        else {
            phi <- plogis(Xphi %*% pars[(nLP+1):(nLP+nPP)] + Xphi.offset)
            phi <- matrix(phi, M, T, byrow=TRUE)
            }
        sigma <- exp(Xdet %*% pars[(nLP+nPP+1):(nLP+nPP+nDP)]+Xdet.offset)
        sigma <- matrix(sigma, M, T, byrow=TRUE)

        switch(mixture,
            P = f <- sapply(k, function(x) dpois(x, lambda)),
            NB = f <- sapply(k, function(x) dnbinom(x, mu=lambda,
                size=exp(pars[nP]))))
        for(i in 1:M) {
            mn <- matrix(0, lk, T)
            for(t in 1:T) {
                if(all(naflag[i,t,]))
                    next
                p <- rep(NA, J)
                switch(survey,
                line = {
                    f.0 <- 2 * dnorm(0, 0, sd=sigma[i, t])
                    int <- 2 * (pnorm(db[-1], 0, sd=sigma[i, t]) -
                        pnorm(db[-(J+1)], 0, sd=sigma[i, t]))
                    p <- int / f.0 / w
                    },
                point = {
                    for(j in 1:J) {
#                        int <- integrate(grhn, db[j], db[j+1],
#                            sigma=sigma[i, t], rel.tol=rel.tol,
#                            stop.on.error=FALSE, subdivisions=50)
#                        if(!identical(int$message, "OK"))
#                            int$value <- NA
                        int <- sigma[i,t]^2 *
                            (1-exp(-db[j+1]^2 / (2*sigma[i,t]^2))) -
                                sigma[i,t]^2 *
                                    (1-exp(-db[j]^2 / (2*sigma[i,t]^2)))
#                        p[j] <- int$value * 2 * pi / a[i,j]
                        p[j] <- int * 2 * pi / a[i,j]
                        }
                    })
                cp <- p * u[i,] * phi[i, t]
                cp[J+1] <- 1 - sum(cp)

                mn[, t] <- lfac.k - lfac.kmyt[i, t,] +
                    sum(y[i, t, !naflag[i,t,]] *
                    log(cp[which(!naflag[i,t,])])) +
                    kmyt[i, t,] * log(cp[J+1])
            }
            g[i,] <- exp(rowSums(mn))
        }
        f[!fin] <- g[!fin] <- 0
        ll <- rowSums(f*g)
        -sum(log(ll))
        }
    },
exp = {
    altdetParms <- paste("rate", colnames(Xdet), sep="")
    if(missing(starts)) {
        starts <- rep(0, nP)
        starts[nLP+nPP+1] <- log(max(db))
        }

    nll <- function(pars) {
        lambda <- exp(Xlam %*% pars[1:nLP] + Xlam.offset)
        if(identical(output, "density"))
            lambda <- lambda * A
        if(T==1)
            phi <- matrix(1, M, T)
        else {
            phi <- plogis(Xphi %*% pars[(nLP+1):(nLP+nPP)] + Xphi.offset)
            phi <- matrix(phi, M, T, byrow=TRUE)
            }
        rate <- exp(Xdet %*% pars[(nLP+nPP+1):(nLP+nPP+nDP)] + Xdet.offset)
        rate <- matrix(rate, M, T, byrow=TRUE)

        switch(mixture,
            P = f <- sapply(k, function(x) dpois(x, lambda)),
            NB = f <- sapply(k, function(x) dnbinom(x, mu=lambda,
                size=exp(pars[nP]))))
        for(i in 1:M) {
            mn <- matrix(0, lk, T)
            for(t in 1:T) {
                if(all(naflag[i,t,]))
                    next
                p <- rep(NA, J)
                switch(survey,
                line = {
                    for(j in 1:J) {
#                        int <- integrate(gxexp, db[j], db[j+1],
#                             rate=rate[i,t], rel.tol=rel.tol,
#                             stop.on.error=FALSE, subdivisions=50)
#                        if(!identical(int$message, "OK"))
#                            int$value <- NA
                        int <- rate[i,t]*(1-exp(-db[j+1]/rate[i,t])) -
                            rate[i,t]*(1-exp(-db[j]/rate[i,t]))
#                        p[j] <- int$value / w[j]
                        p[j] <- int / w[j]
                        }
                    },
                point = {
                    for(j in 1:J) {
                        int <- integrate(grexp, db[j], db[j+1],
                            rate=rate[i, t], rel.tol=rel.tol,
                            stop.on.error=FALSE, subdivisions=50)
                        if(!identical(int$message, "OK"))
                            int$value <- NA
                        p[j] <- int$value * 2 * pi / a[i,j]
                        }
                    })
                cp <- p * u[i,] * phi[i, t]
                cp[J+1] <- 1 - sum(cp)
                mn[, t] <- lfac.k - lfac.kmyt[i, t,] +
                    sum(y[i, t, !naflag[i,t,]] *
                    log(cp[which(!naflag[i,t,])])) +
                    kmyt[i, t,] * log(cp[J+1])
            }
            g[i,] <- exp(rowSums(mn))
        }
        f[!fin] <- g[!fin] <- 0
        ll <- rowSums(f*g)
        -sum(log(ll))
        }
    },
hazard = {
    altdetParms <- paste("shape", colnames(Xdet), sep="")
    if(missing(starts)) {
        starts <- rep(0, nP)
        }
    nll <- function(pars) {
        lambda <- exp(Xlam %*% pars[1:nLP] + Xlam.offset)
        if(identical(output, "density"))
            lambda <- lambda * A

        if(T==1)
            phi <- matrix(1, M, T)
        else {
            phi <- plogis(Xphi %*% pars[(nLP+1):(nLP+nPP)] + Xphi.offset)
            phi <- matrix(phi, M, T, byrow=TRUE)
            }
        shape <- exp(Xdet %*% pars[(nLP+nPP+1):(nLP+nPP+nDP)]+Xdet.offset)
        shape <- matrix(shape, M, T, byrow=TRUE)

        scale <- exp(pars[nLP+nPP+nDP+1])

        switch(mixture,
            P = f <- sapply(k, function(x) dpois(x, lambda)),
            NB = f <- sapply(k, function(x) dnbinom(x, mu=lambda,
                size=exp(pars[nP]))))
        for(i in 1:M) {
            mn <- matrix(0, lk, T)
            for(t in 1:T) {
                if(all(naflag[i,t,]))
                    next
                p <- rep(NA, J)
                switch(survey,
                line = {
                    for(j in 1:J) {
                        int <- integrate(gxhaz, db[j], db[j+1],
                             shape=shape[i,t], scale=scale,
                             rel.tol=rel.tol,
                             stop.on.error=FALSE, subdivisions=50)
                        if(!identical(int$message, "OK"))
                            int$value <- NA
                        p[j] <- int$value / w[j]
                        }
                    },
                point = {
                    for(j in 1:J) {
                        int <- integrate(grhaz, db[j], db[j+1],
                            shape=shape[i, t], scale=scale,
                            rel.tol=rel.tol,
                            stop.on.error=FALSE, subdivisions=50)
                        if(!identical(int$message, "OK"))
                            int$value <- NA
                        p[j] <- int$value * 2 * pi / a[i,j]
                        }
                    })
                cp <- p * u[i,] * phi[i, t]
                cp[J+1] <- 1 - sum(cp)
                mn[, t] <- lfac.k - lfac.kmyt[i, t,] +
                    sum(y[i, t, !naflag[i,t,]] *
                    log(cp[which(!naflag[i,t,])])) +
                    kmyt[i, t,] * log(cp[J+1])
                }
            g[i,] <- exp(rowSums(mn))
            }
        f[!fin] <- g[!fin] <- 0
        ll <- rowSums(f*g)
        -sum(log(ll))
        }
    },
uniform = {
    if(missing(starts)) {
        starts <- rep(0, nP)
        }
    nll <- function(pars) {
        lambda <- exp(Xlam %*% pars[1:nLP] + Xlam.offset)
        if(identical(output, "density"))
            lambda <- lambda * A
        if(T==1)
            phi <- matrix(1, M, T)
        else {
            phi <- plogis(Xphi %*% pars[(nLP+1):(nLP+nPP)] + Xphi.offset)
            phi <- matrix(phi, M, T, byrow=TRUE)
            }
        p <- 1
        switch(mixture,
            P = f <- sapply(k, function(x) dpois(x, lambda)),
            NB = f <- sapply(k, function(x) dnbinom(x, mu=lambda,
                size=exp(pars[nP]))))
        for(i in 1:M) {
            mn <- matrix(0, lk, T)
            for(t in 1:T) {
                if(all(naflag[i,t,]))
                    next
                cp <- p * u[i,] * phi[i, t]
                cp[J+1] <- 1 - sum(cp)
                mn[, t] <- lfac.k - lfac.kmyt[i, t,] +
                    sum(y[i, t, !naflag[i,t,]] *
                    log(cp[which(!naflag[i,t,])])) +
                    kmyt[i, t,] * log(cp[J+1])
            }
            g[i,] <- exp(rowSums(mn))
        }
        f[!fin] <- g[!fin] <- 0
        ll <- rowSums(f*g)
        -sum(log(ll))
    }
})

fm <- optim(starts, nll, method = method, hessian = se, ...)
opt <- fm
if(se) {
  covMat <- tryCatch(solve(fm$hessian), error=function(x)
        simpleError("Hessian is singular. Try using fewer covariates and supplying starting values."))
    if(identical(class(covMat)[1], "simpleError")) {
        warning(covMat$message)
        covMat <- matrix(NA, nP, nP)
        }
    }
else
    covMat <- matrix(NA, nP, nP)
ests <- fm$par
fmAIC <- 2 * fm$value + 2 * nP

names(ests) <- c(lamPars, phiPars, detPars, scalePar, nbPar)

lamEstimates <- unmarkedEstimate(name = "Abundance", short.name = "lambda",
    estimates = ests[1:nLP],
    covMat = as.matrix(covMat[1:nLP, 1:nLP]), invlink = "exp",
    invlinkGrad = "exp")
estimateList <- unmarkedEstimateList(list(lambda=lamEstimates))

if(T>1)
    estimateList@estimates$phi <- unmarkedEstimate(name = "Availability",
        short.name = "phi", estimates = ests[(nLP+1):(nLP+nPP)],
        covMat = as.matrix(covMat[(nLP+1):(nLP+nPP), (nLP+1):(nLP+nPP)]),
        invlink = "logistic", invlinkGrad = "logistic.grad")

if(!identical(keyfun, "uniform"))
    estimateList@estimates$det <- unmarkedEstimate(name = "Detection",
        short.name = "p",
        estimates = ests[(nLP+nPP+1):(nLP+nPP+nDP)],
        covMat = as.matrix(
            covMat[(nLP+nPP+1):(nLP+nPP+nDP), (nLP+nPP+1):(nLP+nPP+nDP)]),
        invlink = "exp", invlinkGrad = "exp")

if(identical(keyfun, "hazard"))
    estimateList@estimates$scale <- unmarkedEstimate(
        name = "Hazard-rate(scale)", short.name = "scale",
        estimates = ests[(nLP+nPP+nDP+1):(nLP+nPP+nDP+1)],
        covMat = covMat[(nLP+nPP+nDP+1):(nLP+nPP+nDP+1),
                        (nLP+nPP+nDP+1):(nLP+nPP+nDP+1), drop=FALSE],
        invlink = "exp", invlinkGrad = "exp")

if(identical(mixture, "NB"))
    estimateList@estimates$alpha <- unmarkedEstimate(name = "Dispersion",
        short.name = "alpha", estimates = ests[nP],
        covMat = as.matrix(covMat[nP, nP]), invlink = "exp",
        invlinkGrad = "exp")

umfit <- new("unmarkedFitGDS", fitType = "gdistsamp",
    call = match.call(), formula = form, formlist = formlist,
    data = data, estimates = estimateList, sitesRemoved = D$removed.sites,
    AIC = fmAIC, opt = opt, negLogLike = fm$value, nllFun = nll,
    mixture=mixture, K=K, keyfun=keyfun, unitsOut=unitsOut, output=output)

return(umfit)
}




