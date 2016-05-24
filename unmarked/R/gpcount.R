
# data will need to be an unmarkedMultFrame
gpcount <- function(lambdaformula, phiformula, pformula, data,
    mixture=c('P', 'NB'), K, starts, method = "BFGS", se = TRUE,
    engine=c('C', 'R'), ...)
{
if(!is(data, "unmarkedFrameGPC"))
    stop("Data is not of class unmarkedFrameGPC.")
mixture <- match.arg(mixture)
engine <- match.arg(engine)
if(identical(mixture, "ZIP") & identical(engine, "R"))
    stop("ZIP mixture not available when 'engine=R'")

formlist <- list(lambdaformula = lambdaformula, phiformula = phiformula,
    pformula = pformula)
form <- as.formula(paste(unlist(formlist), collapse=" "))
D <- getDesign(data, formula = form)

Xlam <- D$Xlam
Xphi <- D$Xphi
Xdet <- D$Xdet
ym <- D$y  # MxJT

Xlam.offset <- D$Xlam.offset
Xphi.offset <- D$Xphi.offset
Xdet.offset <- D$Xdet.offset
if(is.null(Xlam.offset)) Xlam.offset <- rep(0, nrow(Xlam))
if(is.null(Xphi.offset)) Xphi.offset <- rep(0, nrow(Xphi))
if(is.null(Xdet.offset)) Xdet.offset <- rep(0, nrow(Xdet))

if(missing(K) || is.null(K)) {
    K <- max(ym, na.rm=TRUE) + 100
    warning("K was not specified, so was set to max(y)+100 =", K)
}
M <- N <- 0:K
lM <- length(M)
I <- nrow(ym)
T <- data@numPrimary
if(T==1)
    stop("use pcount instead")
J <- numY(data) / T

y <- array(ym, c(I, J, T))

lamPars <- colnames(Xlam)
detPars <- colnames(Xdet)
nLP <- ncol(Xlam)
nPP <- ncol(Xphi)
phiPars <- colnames(Xphi)
nDP <- ncol(Xdet)
nP <- nLP + nPP + nDP + (mixture=='NB')
if(!missing(starts) && length(starts) != nP)
    stop("There should be", nP, "starting values, not", length(starts))

if(identical(engine, "R")) {
# Minus negative log-likelihood
nll <- function(pars) {
    lam <- exp(Xlam %*% pars[1:nLP] + Xlam.offset)
    phi <- plogis(Xphi %*% pars[(nLP+1):(nLP+nPP)] + Xphi.offset)
    phi <- matrix(phi, I, T, byrow=TRUE)
    p <- plogis(Xdet %*% pars[(nLP+nPP+1):(nLP+nPP+nDP)] + Xdet.offset)
    p <- matrix(p, I, byrow=TRUE)
    p <- array(p, c(I, J, T))  # byrow?
    L <- rep(NA, I)
    for(i in 1:I) {
        f <- switch(mixture,
            P = dpois(M, lam[i], log=TRUE),
            NB = dnbinom(M, mu=lam[i], size=exp(pars[nP]), log=TRUE))
#            ZIP = dzip())
        ghi <- rep(0, lM)
        for(t in 1:T) {
            gh <- matrix(-Inf, lM, lM)
            for(m in M) {
                if(m < max(y[i,,], na.rm=TRUE)) {
                    gh[,m+1] <- -Inf
                    next
                }
                if(is.na(phi[i,t])) {
                    g <- rep(0, lM)
                    g[N>m] <- -Inf
                }
                else
                    g <- dbinom(N, m, phi[i,t], log=TRUE)
                h <- rep(0, lM)
                for(j in 1:J) {
                    if(is.na(y[i,j,t]))
                        next
                    h <- h + dbinom(y[i,j,t], N, p[i,j,t], log=TRUE)
                }
                gh[,m+1] <- g + h
            }
            ghi <- ghi + log(colSums(exp(gh))) # sum over N(t)
        }
        fgh <- f + ghi
        L[i] <- sum(exp(fgh)) # sum over M
    }
    return(-sum(log(L)))
}
} else
if(identical(engine, "C")) {
    nll <- function(pars) {
        beta.lam <- pars[1:nLP]
        beta.phi <- pars[(nLP+1):(nLP+nPP)]
        beta.p <- pars[(nLP+nPP+1):(nLP+nPP+nDP)]
        log.alpha <- 1
        if(mixture %in% c("NB", "ZIP"))
            log.alpha <- pars[nP]
        .Call("nll_gpcount",
              ym, Xlam, Xphi, Xdet,
              beta.lam, beta.phi, beta.p, log.alpha,
              Xlam.offset, Xphi.offset, Xdet.offset,
              as.integer(K),
              mixture, T,
              PACKAGE = "unmarked")
    }
}

if(missing(starts)) starts <- rep(0, nP)
fm <- optim(starts, nll, method = method, hessian = se, ...)
opt <- fm
if(se) {
    covMat <- tryCatch(solve(fm$hessian), error=function(x)
        simpleError("Hessian is singular. Try using fewer covariates."))
    if(identical(class(covMat)[1], "simpleError")) {
        warning(covMat$message)
        covMat <- matrix(NA, nP, nP)
        }
    } else covMat <- matrix(NA, nP, nP)
ests <- fm$par
fmAIC <- 2 * fm$value + 2 * nP

if(identical(mixture,"NB")) nbParm <- "alpha"
	else nbParm <- character(0)

names(ests) <- c(lamPars, phiPars, detPars, nbParm)

lamEstimates <- unmarkedEstimate(name = "Abundance", short.name = "lambda",
    estimates = ests[1:nLP],
    covMat = as.matrix(covMat[1:nLP, 1:nLP]), invlink = "exp",
    invlinkGrad = "exp")

phiEstimates <- unmarkedEstimate(name = "Availability",
                                 short.name = "phi",
                                 estimates = ests[(nLP+1):(nLP+nPP)],
                                 covMat = as.matrix(covMat[(nLP+1) :
                                 (nLP+nPP), (nLP+1):(nLP+nPP)]),
                                 invlink = "logistic",
                                 invlinkGrad = "logistic.grad")

detEstimates <- unmarkedEstimate(name = "Detection", short.name = "p",
    estimates = ests[(nLP+nPP+1):(nLP+nPP+nDP)],
    covMat = as.matrix(
        covMat[(nLP+nPP+1):(nLP+nPP+nDP), (nLP+nPP+1):(nLP+nPP+nDP)]),
    invlink = "logistic", invlinkGrad = "logistic.grad")

estimateList <- unmarkedEstimateList(list(lambda=lamEstimates, phi=phiEstimates, det=detEstimates))

if(identical(mixture,"NB"))
    estimateList@estimates$alpha <- unmarkedEstimate(name = "Dispersion",
        short.name = "alpha", estimates = ests[nP],
        covMat = as.matrix(covMat[nP, nP]), invlink = "exp",
        invlinkGrad = "exp")

umfit <- new("unmarkedFitGPC", fitType = "gpcount",
    call = match.call(), formula = form, formlist = formlist,
    data = data, estimates = estimateList, sitesRemoved = D$removed.sites,
    AIC = fmAIC, opt = opt, negLogLike = fm$value, nllFun = nll,
    mixture=mixture, K=K)

return(umfit)
}




