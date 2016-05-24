#### Mallows quasi-likelihood estimator of E. Cantoni and E. Ronchetti (2001)
#### based originally on Eva Cantoni's S-plus code "robGLM"

## FIXME{MM}: All these expression()s  and  eval()s -- once were really slick and fast.
## -----  Nowadays, with 'codetools' and the byte-compiler, they "just don't fit anymore"
## including those globalVariables() {also in other places!}:
globalVariables(c("residP", "residPS", "dmu.deta", "snu"), add=TRUE)

##' @title
##' @param wts a character string \dQuote{weights.on.x} specifying how weights should be computed
##'            *or* a numeric vector of final weights in which case nothing is computed.
##' @param X  n x p  design matrix aka model.matrix()
##' @param intercept logical, if true, X[,] has an intercept column which should
##'                  not be used for rob.wts
##' @return n-vector of non-negative weights
##' @author Martin Maechler
robXweights <- function(wts, X, intercept=TRUE) {
    stopifnot(length(d <- dim(X)) == 2, is.logical(intercept))
    nobs <- d[1]
    if(d[2]) { ## X has >= 1 column, and hence there *are* coefficients in the end
        if(is.character(wts)){
	    switch(wts,
		   "none" = rep.int(1, nobs),
		   "hat" = wts_HiiDist(X)^2, # = (1 - Hii)^2
		   "robCov" = wts_RobDist(X, intercept, covFun = MASS::cov.rob),
		   ## MCD is currently problematic: many singular subsamples
		   "covMcd" = wts_RobDist(X, intercept, covFun = covMcd),
		   stop("Weighting method", sQuote(wts),
			" is not implemented"))
	}
	## (new; 2013-07-05; -> robustbase 0.9-9)
	else if(is.list(wts)) {
	    if(length(wts) == 1 && is.function(covF <- wts[[1]]))
		wts_RobDist(X, intercept, covFun = covF)
	    else stop("if a list, weights.on.x must contain a covariance function such as covMcd()")
	}
	else if(is.function(wts)) {
	    wts(X, intercept)
	}
	else {
	    if(!is.numeric(wts) || length(wts) != nobs)
		## FIXME: "when not a string, a list, or a function, then ..."
		stop(gettextf("weights.on.x needs %d none-negative values",
			      nobs), domain=NA)
            if(any(wts) < 0)
                stop("All weights.on.x must be none negative")
        }
    }
    else ## p = ncoef == 0 {maybe intercept, but that's not relevant here}
        rep.int(1,nobs)
}


##' @param intercept logical, if true, X[,] has an intercept column which should
##'                  not be used for rob.wts
glmrobMqle <-
    function(X, y, weights = NULL, start = NULL, offset = NULL,
	     family, weights.on.x = "none",
	     control = glmrobMqle.control(), intercept = TRUE,
             trace = FALSE)
{
    ## To DO:
    ## o weights are not really implemented as *extra* user weights; rather as "glm-weights"
    ## o offset is not fully implemented (really? -- should have test case!)

    if(!is.matrix(X)) X <- as.matrix(X)
## never used:
##     xnames <- dimnames(X)[[2]]
##     ynames <- if (is.matrix(y)) rownames(y) else names(y)
    nobs <- NROW(y)
    stopifnot(nobs == nrow(X))
    if (is.null(weights))
	weights <- rep.int(1, nobs)
    else if(any(weights <= 0))
	stop("All weights must be positive")
    if (is.null(offset))
	offset <- rep.int(0, nobs) else if(!all(offset==0))
	    warning("'offset' not fully implemented")
    variance <- family$variance
    linkinv <- family$linkinv
    if (!is.function(variance) || !is.function(linkinv))
	stop("illegal 'family' argument")
    mu.eta <- family$mu.eta
    if (is.null(valideta <- family$valideta)) valideta <- function(eta) TRUE
    if (is.null(validmu	 <- family$validmu))  validmu <-  function(mu) TRUE

    ncoef <- ncol(X)
    w.x <- robXweights(weights.on.x, X=X, intercept=intercept)

### Initializations
    stopifnot(control$maxit >= 1, (tcc <- control$tcc) >= 0)

    ## note that etastart and mustart are used to make 'family$initialize' run
    etastart <- NULL;  mustart <- NULL
    ## note that 'weights' are used and set by binomial()$initialize !
    eval(family$initialize) ## --> n, mustart, y and weights (=ni)
    ni <- as.vector(weights)# dropping attributes for computation
    ##
    if(is.null(start))
	start <- glm.fit(x = X, y = y, weights = weights, offset = offset,
			 family = family)$coefficients
    if(any(ina <- is.na(start))) {
	cat("initial start 'theta' has NA's; eliminating columns X[, j];",
	    "j = ", pasteK(which(ina)),"\n")
	theta.na <- start
	X <- X[, !ina, drop = FALSE]
	start <- glm.fit(x = X, y = y, weights = weights, offset = offset,
			 family = family)$coefficients
	if(any(is.na(start)))
	    stop("start 'theta' has still NA's .. badly singular x\n")
	## FIXME
	ncoef <- length(start)
    }

    thetaOld <- theta <- as.vector(start) # as.v*(): dropping attributes
    eta <- as.vector(X %*% theta)
    mu <- linkinv(eta) # mu estimates pi (in [0,1]) at the binomial model
    if (!(validmu(mu) && valideta(eta)))
	stop("Cannot find valid starting values: You need help")
    ##
    switch(family$family,
	   "binomial" = {
	       Epsi.init <- EpsiBin.init
	       Epsi <- EpsiBin
	       EpsiS <- EpsiSBin
	       Epsi2 <- Epsi2Bin
               phiEst <- phiEst.cl <- 1
	   },
	   "poisson" = {
	       Epsi.init <- EpsiPois.init
	       Epsi <- EpsiPois
	       EpsiS <- EpsiSPois
	       Epsi2 <- Epsi2Pois
               phiEst <- phiEst.cl <- expression({1})
	   },
           "gaussian" = {
               Epsi.init <- EpsiGaussian.init
               Epsi <- EpsiGaussian
               EpsiS <- EpsiSGaussian
               Epsi2 <- Epsi2Gaussian
               phiEst.cl <- phiGaussianEst.cl
               phiEst <- phiGaussianEst
           },
          "Gamma" = { ## added by ARu
             Epsi.init <- EpsiGamma.init
             Epsi <- EpsiGamma
             EpsiS <- EpsiSGamma
             Epsi2 <- Epsi2Gamma
             phiEst.cl <- phiGammaEst.cl
             phiEst <- phiGammaEst
           },
           ## else
           stop(gettextf("family '%s' not yet implemented", family$family),
                domain=NA)
	   )

    sV <- NULL # FIXME workaround for codetools

    comp.V.resid <- expression({
	Vmu <- variance(mu)
	if (any(is.na(Vmu)))  stop("NAs in V(mu)")
	if (any(Vmu == 0))    stop("0s in V(mu)")
	sVF <- sqrt(Vmu)   # square root of variance function
	residP <- (y - mu)* sni/sVF  # Pearson residuals
    })

    comp.scaling <- expression({
      sV <- sVF * sqrt(phi)
      residPS <- residP/sqrt(phi) # scaled Pearson residuals
    })

    comp.Epsi.init <- expression({
	## d mu / d eta :
	dmu.deta <- mu.eta(eta)
	if (any(is.na(dmu.deta))) stop("NAs in d(mu)/d(eta)")
	## "Epsi init" :
	H <- floor(mu*ni - tcc* sni*sV)
	K <- floor(mu*ni + tcc* sni*sV)
	eval(Epsi.init)
    })


### Iterations

    if(trace && ncoef) {
        cat("Initial theta: \n")
        local({names(theta) <- names(start); print(theta) })

        digits <- max(1, getOption("digits") - 5)
	w.th.1 <- 6+digits # width of one number; need 8 for 2 digits: "-4.8e-11"
	width.th <- ncoef*(w.th.1 + 1) - 1
	cat(sprintf("%3s | %*s | %12s\n",
		    "it", width.th, "d{theta}", "rel.change"))
	mFormat <- function(x, wid) {
	    r <- formatC(x, digits=digits, width=wid)
	    sprintf("%*s", wid, sub("e([+-])0","e\\1", r))
	}
    }

    sni <- sqrt(ni)
    eval(comp.V.resid) #-> (Vmu, sVF, residP)
    phi <- eval(phiEst.cl)
    ## Determine the range of phi values based on the distribution of |residP|
    Rphi <- c(1e-12, 3*median(abs(residP)))^2
    conv <- FALSE
    if(ncoef) for (nit in 1:control$maxit) {
        eval(comp.scaling) #-> (sV, residPS)
        eval(comp.Epsi.init)
	## Computation of alpha and (7) using matrix column means:
	cpsi <- pmax.int(-tcc, pmin.int(residPS,tcc)) - eval(Epsi)
	EEq <- colMeans(cpsi * w.x * sni/sV * dmu.deta * X)
	##
	## Solve  1/n (t(X) %*% B %*% X) %*% delta.coef	  = EEq
	DiagB <- eval(EpsiS) /(sni*sV) * w.x * (ni*dmu.deta)^2
        if(any(n0 <- ni == 0)) DiagB[n0] <- 0 # instead of NaN
	Dtheta <- solve(crossprod(X, DiagB*X)/nobs, EEq)
	if (any(!is.finite(Dtheta))) {
	    warning("Non-finite coefficients at iteration ", nit)
	    break
	}
	theta <- thetaOld + Dtheta
	eta <- as.vector(X %*% theta) + offset
	mu <- linkinv(eta)

        ## estimation of the dispersion parameter
        eval(comp.V.resid)
        phi <- eval(phiEst)

	## Check convergence: relative error < tolerance
	relE <- sqrt(sum(Dtheta^2)/max(1e-20, sum(thetaOld^2)))
	conv <- relE <= control$acc
        if(trace) {
            cat(sprintf("%3d | %*s | %12g\n", nit, width.th,
                        paste(mFormat(Dtheta, w.th.1),
                              collapse=" "), relE))
        }
	if(conv)
	    break
	thetaOld <- theta
    } ## end of iteration
    else { ## ncoef == 0
	conv <- TRUE
	nit <- 0
    }
    if (!conv)
	warning("Algorithm did not converge")

    eps <- 10 * .Machine$double.eps
    switch(family$family,
	   "binomial" = {
	       if (any(mu/weights > 1 - eps) || any(mu/weights < eps))
		   warning("fitted probabilities numerically 0 or 1 occurred")
	   },
	   "poisson" = {
	       if (any(mu < eps))
		   warning("fitted rates numerically 0 occurred")
	   })

    eval(comp.V.resid) #-> (Vmu, sVF, residP)
    eval(comp.scaling) #-> (sV, residPS)

    ## Estimated asymptotic covariance of the robust estimator
    if(ncoef) {
	eval(comp.Epsi.init)
	alpha <- colMeans(eval(Epsi) * w.x * sni/sV * dmu.deta * X)
	DiagA <- eval(Epsi2) / (ni*sV^2)* w.x^2* (ni*dmu.deta)^2
	matQ  <- crossprod(X, DiagA*X)/nobs - tcrossprod(alpha, alpha)

	DiagB <- eval(EpsiS) / (sni*sV)* w.x * (ni*dmu.deta)^2
        if(any(n0 <- ni == 0)) DiagB[n0] <- 0 # instead of NaN
	matM <- crossprod(X, DiagB*X)/nobs
	matMinv <- solve(matM)
	asCov <-  matMinv %*% matQ %*% matMinv / nobs
    } else { ## ncoef == 0
	matM <- matQ <- asCov <- matrix(, 0,0)
    }

    if(any(ina)) {# put NA's back, extending theta[] to "original length"
	ok <- !ina
	theta.na[ok] <- theta ; theta <- theta.na
	## also extend the "p x p" matrices with NA's --
	##No : lm() and glm() also do *not* do this
	##No  p <- length(theta)
	##No  nm <- names(theta)
	##No  M <- matrix(as.numeric(NA), p, p, dimnames = list(nm,nm))
	##No  Mn <- M; Mn[ok, ok] <- asCov ; asCov <- Mn
	##No  Mn <- M; Mn[ok, ok] <- matM  ; matM  <- Mn
	##No  Mn <- M; Mn[ok, ok] <- matQ  ; matQ  <- Mn
    }

    w.r <- pmin(1, tcc/abs(residPS))
    names(mu) <- names(eta) <- names(residPS) # re-add after computation
    list(coefficients = theta, residuals = residP, # s.resid = residPS,
         fitted.values = mu,
	 w.r = w.r, w.x = w.x, ni = ni, dispersion = phi, cov = asCov,
         matM = matM, matQ = matQ, tcc = tcc, family = family,
         linear.predictors = eta, deviance = NULL, iter = nit, y = y,
         converged = conv)
}


## NB: X  is model.matrix() aka design matrix used; typically including an intercept
wts_HiiDist <- function(X) {
    ## Hii := diag( tcrossprod( qr.Q(qr(X)) ) ) == rowSums( qr.Q(qr(X)) ^2 ) :
    x <- qr(X)
    Hii <- rowSums(qr.qy(x, diag(1, nrow = NROW(X), ncol = x$rank))^2)
    (1-Hii)
}

##' Compute robustness weights depending on the design 'X' only,
##' using robust(ified) Mahalanobis distances.
##' This is an auxiliary function for robXweights() activated typically by
##' weights.on.x = "..." from regression functions
##' @title Compute Robust Weights based on Robustified Mahalanobis - Distances
##' @param X n x p  numeric matrix
##' @param intercept logical; should be true iff  X[,1] is a column with the intercept
##' @param covFun function for computing a \bold{robust} covariance matrix;
##'        e.g., MASS::cov.rob(), or covMcd().
##' @return n-vector of non-negative weights.
##' @author Martin Maechler
wts_RobDist <- function(X, intercept, covFun)
{
    D2 <- if(intercept) { ## X[,] has intercept column which should not be used for rob.wts
	X <- X[, -1, drop=FALSE]
	Xrc <- covFun(X)
	mahalanobis(X, center = Xrc$center, cov = Xrc$cov)
    }
    else { ## X[,]  can be used directly
	if(!is.matrix(X)) X <- as.matrix(X)
	Xrc <- covFun(X)
	S <- Xrc$cov + tcrossprod(Xrc$center)
	mahalanobis(X, center = FALSE, cov = S)
    }
    p <- ncol(X) ## E[chi^2_p] = p
    1/sqrt(1+ pmax.int(0, 8*(D2 - p)/sqrt(2*p)))
}


## MM: 'acc' seems a misnomer to me, but it's inherited from  MASS::rlm
glmrobMqle.control <-
    function(acc = 1e-04, test.acc = "coef", maxit = 50, tcc = 1.345)
{
    if (!is.numeric(acc) || acc <= 0)
	stop("value of acc must be > 0")
    if (test.acc != "coef")
	stop("Only 'test.acc = \"coef\"' is currently implemented")
    ## if (!(any(test.vec == c("coef", "resid"))))
    ##	  stop("invalid argument for test.acc")
    if (!is.numeric(maxit) || maxit <= 0)
	stop("maximum number of iterations must be > 0")
    if (!is.numeric(tcc) || tcc <= 0)
	stop("value of the tuning constant c (tcc) must be > 0")
    list(acc = acc, test.acc = test.acc, maxit = maxit, tcc = tcc)
}


### ----------------- E[ f(psi ( X ) ) ] -------------------------------

## MM: These are now expressions instead of functions
##   since 'Epsi*' and 'Epsi2*' are *always* called together
##   and 'EpsiS*' when called is always after the other two
## ==> do common computations only once in Epsi*.init  ==> more efficient!
##
##   FIXME(2): Some of these fail when Huber's "c", 'tcc' is = +Inf
##   -----    --> ../../robGLM1/R/rglm.R

## FIXME:  Do use a "robFamily", a  *list* of functions
## ------  which all have the same environment
##   ===> can get same efficiency as expressions, but better OOP


### --- Poisson -- family ---

EpsiPois.init <- expression(
{
    dpH <- dpois(H, mu); dpH1 <- dpois(H-1, mu)
    dpK <- dpois(K, mu); dpK1 <- dpois(K-1, mu)
    pHm1 <- ppois(H-1, mu) ; pH <- pHm1 + dpH # = ppois(H,*)
    pKm1 <- ppois(K-1, mu) ; pK <- pKm1 + dpK # = ppois(K,*)
    E2f <- mu*(dpH1 - dpH - dpK1 + dpK) + pKm1 - pHm1
})

EpsiPois <- expression(
{
    tcc*(1 - pK - pH) + mu*(dpH - dpK)/sV
})

Epsi2Pois <- expression(
{
    ## Calculation of E(psi^2) for the diagonal elements of A in matrix Q:
    tcc^2 * (pH + 1 - pK) + E2f
})

EpsiSPois <- expression(
{
    ## Calculation of E(psi*s) for the diagonal elements of B in the
    ## expression matrix M = 1/n t(X) %*% B %*% X:
    tcc*(dpH + dpK) + E2f / sV
})


### --- Binomial -- family ---

EpsiBin.init <- expression({
    pK <- pbinom(K, ni, mu)
    pH <- pbinom(H, ni, mu)
    pKm1 <- pbinom(K-1, pmax.int(0, ni-1), mu)
    pHm1 <- pbinom(H-1, pmax.int(0, ni-1), mu)
    pKm2 <- pbinom(K-2, pmax.int(0, ni-2), mu)
    pHm2 <- pbinom(H-2, pmax.int(0, ni-2), mu)

    ## QlV = Q / V, where Q = Sum_j (j - mu_i)^2 * P[Y_i = j]
    ## i.e.  Q =	     Sum_j j(j-1)* P[.] +
    ##		 (1- 2*mu_i) Sum_j   j	 * P[.] +
    ##		     mu_i^2  Sum_j	   P[.]
    QlV <- mu/Vmu*(mu*ni*(pK-pH) +
		   (1 - 2*mu*ni) * ifelse(ni == 1, (H <= 0)*(K >= 1), pKm1 - pHm1) +
		   (ni - 1) * mu * ifelse(ni == 2, (H <= 1)*(K >= 2), pKm2 - pHm2))
})

EpsiBin <- expression(
{
    tcc*(1 - pK - pH) +
	ifelse(ni == 1, (- (H < 0) + (K >= 1) ) * sV,
	       (pKm1 - pHm1 - pK + pH) * mu * sni/sV)
})

Epsi2Bin <- expression(
{
    ## Calculation of E(psi^2) for the diagonal elements of A in matrix Q:
    tcc^2*(pH + 1 - pK) + QlV
})

EpsiSBin <- expression(
{
    ## Calculation of E(psi*s) for the diagonal elements of B in the
    ## expression matrix M = (X' B X)/n
    mu/Vmu*(tcc*(pH - ifelse(ni == 1, H >= 1, pHm1)) +
	    tcc*(pK - ifelse(ni == 1, K > 0,  pKm1))) + ifelse(ni == 0, 0, QlV / (sni*sV))
})

### --- Gaussian -- family ---

EpsiGaussian.init <- expression({
    dc <- dnorm(tcc)
    pc <- pnorm(tcc)
})

EpsiGaussian <- expression( 0 )

EpsiSGaussian <- expression( 2*pc-1 )

Epsi2Gaussian <- expression( 2*tcc^2*(1-pc)-2*tcc*dc+2*pc-1 )

phiGaussianEst.cl <- expression(
{
    ## Classical estimation of the dispersion paramter phi = sigma^2
    sum(((y - mu)/mu)^2)/(nobs - ncoef)
})

phiGaussianEst <- expression(
{
    sphi <- mad(residP, center=0)^2
})

### --- Gamma -- family ---

Gmn <- function(t, nu) {
    ## Gm corrresponds to G * nu^((nu-1)/2) / Gamma(nu)
    snu <- sqrt(nu)
    snut <- snu+t
    r <- numeric(length(snut))
    ok <- snut > 0
    r[ok] <- {
	nu <- nu[ok]; snu <- snu[ok]; snut <- snut[ok]
	exp((nu-1)/2*log(nu) - lgamma(nu) - snu*snut + nu*log(snut))
    }
    r
}

EpsiGamma.init <- expression({

    nu <- 1/phi      ## form parameter nu
    snu <- 1/sqrt(phi) ## == sqrt (nu)

    pPtc <- pgamma(snu + c(-tcc,tcc), shape=nu, rate=snu)
    pMtc <- pPtc[1]
    pPtc <- pPtc[2]

    aux2 <- tcc*snu
    GLtcc <- Gmn(-tcc,nu)
    GUtcc <- Gmn( tcc,nu)
})

EpsiGamma <- expression( tcc*(1-pPtc-pMtc) + GLtcc - GUtcc )

EpsiSGamma <- expression( ((GLtcc - GUtcc) + snu*(pPtc-pMtc))/mu )

Epsi2Gamma <- expression({
    (tcc^2*(pMtc+1-pPtc) + (pPtc-pMtc) +
     (GLtcc*(1-aux2) - GUtcc*(1+aux2))/snu )
})


phiGammaEst.cl <- expression(
{
    ## Classical moment estimation of the dispersion parameter phi
    sum(((y - mu)/mu)^2)/(nobs-ncoef)
})

phiGammaEst <- expression(
{
    ## robust estimation of the dispersion parameter by
    ## Huber's proposal 2
    sphi <- uniroot(Huberprop2, interval=Rphi,
                    ns.resid=residP, mu=mu, Vmu=Vmu, tcc=tcc)$root
})

Huberprop2 <- function(phi, ns.resid, mu, Vmu, tcc)
{
    eval(EpsiGamma.init)
    compEpsi2 <- eval(Epsi2Gamma)
    nobs <- length(mu)
    ## return h :=
    sum(pmax.int(-tcc, pmin.int(ns.resid*snu, tcc))^2) -  nobs*compEpsi2
}

if(FALSE) ## no-eval version
Huberprop2 <- function(phi, ns.resid, mu, Vmu, tcc)
{
    nobs <- length(mu)
    nu <- 1/phi         ## form parameter  nu
    snu <- 1/sqrt(phi)  ## sqrt (nu)
    pPtc <- pgamma(snu + c(-tcc,tcc), shape=nu, rate=snu)
    pMtc <- pPtc[1]
    pPtc <- pPtc[2]

    ts <- tcc*snu
    GLtcc <- Gmn(-tcc,nu) *(1-ts)/snu
    GUtcc <- Gmn( tcc,nu) *(1+ts)/snu
    ##
    compEpsi2 <- tcc^2 + (pPtc - pMtc)*(1-tcc^2) + GLtcc - GUtcc
    ## return h :=
    sum(pmax.int(-tcc, pmin.int(ns.resid*snu, tcc))^2) -  nobs*compEpsi2
}
