


pcountOpen <- function(lambdaformula, gammaformula, omegaformula, pformula,
    data, mixture=c("P", "NB", "ZIP"), K, dynamics=c("constant", "autoreg",
    "notrend", "trend", "ricker", "gompertz"),
    fix=c("none", "gamma", "omega"),
    starts, method="BFGS", se=TRUE, immigration=FALSE, iotaformula=~1, ...)
{
mixture <- match.arg(mixture)
dynamics <- match.arg(dynamics)
fix <- match.arg(fix)
## if(identical(dynamics, "notrend") &
##    !identical(lambdaformula, omegaformula))
##     stop("lambdaformula and omegaformula must be identical for notrend model")
if((identical(dynamics, "constant") || identical(dynamics, "notrend")) & immigration)
    stop("You can not include immigration in the constant or notrend models")
formlist <- list(lambdaformula=lambdaformula, gammaformula=gammaformula,
    omegaformula=omegaformula, pformula=pformula, iotaformula=iotaformula)
formula <- as.formula(paste(unlist(formlist), collapse=" "))
D <- getDesign(data, formula)
y <- D$y

Xlam <- D$Xlam
Xgam <- D$Xgam
Xom <- D$Xom
Xp <- D$Xp
Xiota <- D$Xiota

delta <- D$delta; go.dims <- D$go.dims
deltamax <- max(delta, na.rm=TRUE)
M <- nrow(y)
T <- data@numPrimary
J <- ncol(getY(data)) / T

Xlam.offset <- D$Xlam.offset
Xgam.offset <- D$Xgam.offset
Xom.offset <- D$Xom.offset
Xp.offset <- D$Xp.offset
Xiota.offset <- D$Xiota.offset
if(is.null(Xlam.offset)) Xlam.offset <- rep(0, M)
if(is.null(Xgam.offset)) Xgam.offset <- rep(0, M*(T-1))
if(is.null(Xom.offset)) Xom.offset <- rep(0, M*(T-1))
if(is.null(Xp.offset)) Xp.offset <- rep(0, M*T*J)
if(is.null(Xiota.offset)) Xiota.offset <- rep(0, M*(T-1))

yna <- is.na(y)
yna[] <- as.integer(yna)
y <- array(y, c(M, J, T))
ytna <- apply(is.na(y), c(1,3), all)
ytna <- matrix(ytna, nrow=M)
ytna[] <- as.integer(ytna)

first <- apply(!ytna, 1, function(x) min(which(x)))
last  <- apply(!ytna, 1, function(x) max(which(x)))

if(missing(K)) {
    K <- max(y, na.rm=T) + 20
    warning("K was not specified and was set to ", K, ".")
}
if(K <= max(y, na.rm = TRUE))
    stop("specified K is too small. Try a value larger than any observation")
k <- 0:K
lk <- length(k)

lamParms <- colnames(Xlam)
gamParms <- colnames(Xgam)
omParms <- colnames(Xom)
detParms <- colnames(Xp)
nAP <- ncol(Xlam)
nGP <- ncol(Xgam)
nOP <- ncol(Xom)
nDP <- ncol(Xp)

if(immigration) {
  iotaParms <- colnames(Xiota)
  nIP <- ncol(Xiota)
} else {
  nIP <- 0
  iotaParms <- character(0)
}
if(identical(fix, "gamma")) {
    if(!identical(dynamics, "constant"))
        stop("dynamics must be constant when fixing gamma or omega")
    if(nGP > 1)
        stop("gamma covariates not allowed when fix==gamma")
    else {
        nGP <- 0
        gamParms <- character(0)
    }
}
else if(identical(dynamics, "notrend")) {
    if(nGP > 1)
        stop("gamma covariates not allowed when dyamics==notrend")
    else {
        nGP <- 0
        gamParms <- character(0)
    }
}

if(identical(fix, "omega")) {
    if(!identical(dynamics, "constant"))
        stop("dynamics must be constant when fixing gamma or omega")
    if(nOP > 1)
        stop("omega covariates not allowed when fix==omega")
    else {
        nOP <- 0
        omParms <- character(0)
    }
} else if(identical(dynamics, "trend")) {
    if(nOP > 1)
        stop("omega covariates not allowed when dynamics='trend'")
    else {
        nOP <- 0
        omParms <- character(0)
    }
}


nP <- nAP + nGP + nOP + nDP + nIP + (mixture!="P")
if(!missing(starts) && length(starts) != nP)
    stop(paste("The number of starting values should be", nP))

ym <- matrix(y, nrow=M)


# Create indices (should be written in C++)
I <- cbind(rep(k, times=lk),
           rep(k, each=lk))
I1 <- I[I[,1] <= I[,2],]
Ib <- Ip <- list()
for(i in 1:nrow(I)) {
    Z <- 0:min(I[i,])
    Ib[[i]] <- which((I1[,1] %in% Z) & (I1[,2]==I[i,1])) - 1
    Ip[[i]] <- as.integer(I[i,2]-Z)
}




nll <- function(parms) {
    beta.lam <- parms[1:nAP]
    beta.gam <- parms[(nAP+1):(nAP+nGP)]
    beta.om <- parms[(nAP+nGP+1):(nAP+nGP+nOP)]
    beta.p <- parms[(nAP+nGP+nOP+1):(nAP+nGP+nOP+nDP)]
    beta.iota <- parms[(nAP+nGP+nOP+nDP+1):(nAP+nGP+nOP+nDP+nIP)]
    log.alpha <- 1
    if(mixture %in% c("NB", "ZIP"))
        log.alpha <- parms[nP]
    .Call("nll_pcountOpen",
          ym,
          Xlam, Xgam, Xom, Xp, Xiota,
          beta.lam, beta.gam, beta.om, beta.p, beta.iota, log.alpha,
          Xlam.offset, Xgam.offset, Xom.offset, Xp.offset, Xiota.offset,
          ytna, yna,
          lk, mixture, first, last, M, J, T,
          delta, dynamics, fix, go.dims, immigration,
          I, I1, Ib, Ip,
          PACKAGE = "unmarked")
}
if(missing(starts))
    starts <- rep(0, nP)
#browser()
fm <- optim(starts, nll, method=method, hessian=se, ...)

opt <- fm
ests <- fm$par
if(identical(mixture,"NB"))
    nbParm <- "alpha"
else if(identical(mixture, "ZIP"))
    nbParm <- "psi"
else
    nbParm <- character(0)
names(ests) <- c(lamParms, gamParms, omParms, detParms, iotaParms, nbParm)
if(se) {
	covMat <- tryCatch(solve(fm$hessian), error=function(x)
		simpleError("Hessian is not invertible. Try using fewer covariates or providing starting values."))
	if(class(covMat)[1] == "simpleError") {
		print(covMat$message)
		covMat <- matrix(NA, nP, nP)
		}
    } else covMat <- matrix(NA, nP, nP)

fmAIC <- 2*fm$value + 2*nP

lamEstimates <- unmarkedEstimate(name = "Abundance", short.name = "lam",
    estimates = ests[1:nAP], covMat = as.matrix(covMat[1:nAP,1:nAP]),
    invlink = "exp", invlinkGrad = "exp")
detEstimates <- unmarkedEstimate(name = "Detection", short.name = "p",
    estimates = ests[(nAP+nGP+nOP+1) : (nAP+nGP+nOP+nDP)],
    covMat = as.matrix(covMat[(nAP+nGP+nOP+1) : (nAP+nGP+nOP+nDP),
        (nAP+nGP+nOP+1) : (nAP+nGP+nOP+nDP)]),
        invlink = "logistic", invlinkGrad = "logistic.grad")
estimateList <- unmarkedEstimateList(list(lambda=lamEstimates))
gamName <- switch(dynamics,
                  constant = "gamConst",
                  autoreg = "gamAR",
                  notrend = "",
                  trend = "gamTrend",
                  ricker = "gamRicker",
                  gompertz = "gamGomp")
if(!(identical(fix, "gamma") | identical(dynamics, "notrend")))
    estimateList@estimates$gamma <- unmarkedEstimate(name =
        ifelse(identical(dynamics, "constant") | identical(dynamics, "autoreg"),
        "Recruitment", "Growth Rate"), short.name = gamName,
        estimates = ests[(nAP+1) : (nAP+nGP)], covMat = as.matrix(covMat[(nAP+1) :
                           (nAP+nGP), (nAP+1) : (nAP+nGP)]),
        invlink = "exp", invlinkGrad = "exp")
if(!(identical(fix, "omega") | identical(dynamics, "trend"))) {
  if(identical(dynamics, "constant") | identical(dynamics, "autoreg") | identical(dynamics, "notrend"))
    estimateList@estimates$omega <- unmarkedEstimate(
        name="Apparent Survival",
        short.name = "omega", estimates = ests[(nAP+nGP+1) :(nAP+nGP+nOP)],
        covMat = as.matrix(covMat[(nAP+nGP+1) : (nAP+nGP+nOP),
            (nAP+nGP+1) : (nAP+nGP+nOP)]),
        invlink = "logistic", invlinkGrad = "logistic.grad")
  else if(identical(dynamics, "ricker"))
    estimateList@estimates$omega <- unmarkedEstimate(
        name="Carrying Capacity",
        short.name = "omCarCap", estimates = ests[(nAP+nGP+1) :(nAP+nGP+nOP)],
        covMat = as.matrix(covMat[(nAP+nGP+1) : (nAP+nGP+nOP),
            (nAP+nGP+1) : (nAP+nGP+nOP)]),
        invlink = "exp", invlinkGrad = "exp")
  else
    estimateList@estimates$omega <- unmarkedEstimate(
        name="Carrying Capacity",
        short.name = "omCarCap", estimates = ests[(nAP+nGP+1) :(nAP+nGP+nOP)],
        covMat = as.matrix(covMat[(nAP+nGP+1) : (nAP+nGP+nOP),
            (nAP+nGP+1) : (nAP+nGP+nOP)]),
        invlink = "exp", invlinkGrad = "exp")
  }
estimateList@estimates$det <- detEstimates
if(immigration) {
    estimateList@estimates$iota <- unmarkedEstimate(
        name="Immigration",
        short.name = "iota", estimates = ests[(nAP+nGP+nOP+nDP+1) :(nAP+nGP+nOP+nDP+nIP)],
        covMat = as.matrix(covMat[(nAP+nGP+nOP+nDP+1) : (nAP+nGP+nOP+nDP+nIP),
            (nAP+nGP+nOP+nDP+1) : (nAP+nGP+nOP+nDP+nIP)]),
        invlink = "exp", invlinkGrad = "exp")
}
if(identical(mixture, "NB")) {
    estimateList@estimates$alpha <- unmarkedEstimate(name = "Dispersion",
        short.name = "alpha", estimates = ests[nP],
        covMat = as.matrix(covMat[nP, nP]), invlink = "exp",
        invlinkGrad = "exp")
    }
if(identical(mixture, "ZIP")) {
    estimateList@estimates$psi <- unmarkedEstimate(name = "Zero-inflation",
        short.name = "psi", estimates = ests[nP],
        covMat = as.matrix(covMat[nP, nP]), invlink = "logistic",
        invlinkGrad = "logistic.grad")
    }
umfit <- new("unmarkedFitPCO", fitType = "pcountOpen",
    call = match.call(), formula = formula, formlist = formlist, data = data,
    sitesRemoved=D$removed.sites, estimates = estimateList, AIC = fmAIC,
    opt = opt, negLogLike = fm$value, nllFun = nll, K = K, mixture = mixture,
    dynamics = dynamics, immigration = immigration)
return(umfit)
}









