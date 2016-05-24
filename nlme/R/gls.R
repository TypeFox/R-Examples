###  Fit a linear model with correlated errors and/or heteroscedasticity
###
### Copyright 1997-2003  Jose C. Pinheiro,
###                      Douglas M. Bates <bates@stat.wisc.edu>
### Copyright 2005-2016  The R Core team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

gls <-
    ## fits linear model with serial correlation and variance functions,
    ## by maximum likelihood using a Newton-Raphson algorithm.
    function(model,
             data = sys.frame(sys.parent()),
             correlation = NULL,
             weights = NULL,
             subset,
             method = c("REML", "ML"),
             na.action = na.fail,
             control = list(),
             verbose = FALSE)
{
    Call <- match.call()
    ## control parameters
    controlvals <- glsControl()
    if (!missing(control)) controlvals[names(control)] <- control

    ##
    ## checking arguments
    ##
    if (!inherits(model, "formula") || length(model) != 3L) {
        stop("\nmodel must be a formula of the form \"resp ~ pred\"")
    }
    method <- match.arg(method)
    REML <- method == "REML"
    ## check if correlation is present and has groups
    groups <- if (!is.null(correlation)) getGroupsFormula(correlation) ## else NULL
    ## create a gls structure containing the plug-ins
    glsSt <-
        glsStruct(corStruct = correlation, varStruct = varFunc(weights))

    ## we need to resolve '.' in the formula here
    model <- terms(model, data=data)
    ## extract a data frame with enough information to evaluate
    ## formula, groups, corStruct, and varStruct
    mfArgs <- list(formula = asOneFormula(formula(glsSt), model, groups),
                   data = data, na.action = na.action)
    if (!missing(subset)) {
        mfArgs[["subset"]] <- asOneSidedFormula(Call[["subset"]])[[2L]]
    }
    mfArgs$drop.unused.levels <- TRUE
    dataMod <- do.call(model.frame, mfArgs)
    origOrder <- row.names(dataMod)	# preserve the original order
    if (!is.null(groups)) {
        ## sort the model.frame by groups and get the matrices and parameters
        ## used in the estimation procedures
        ## always use innermost level of grouping
        groups <- eval(substitute(~ 1 | GR, list(GR = groups[[2L]])))
        grps <- getGroups(dataMod, groups,
                          level = length(getGroupsFormula(groups, asList = TRUE)))
        ## ordering data by groups
        ord <- order(grps)
        grps <- grps[ord]
        dataMod <- dataMod[ord, ,drop = FALSE]
        revOrder <- match(origOrder, row.names(dataMod)) # putting in orig. order
    } else grps <- NULL

    ## obtaining basic model matrices
    X <- model.frame(model, dataMod)
    ## keeping the contrasts for later use in predict
    contr <- lapply(X, function(el)
                    if (inherits(el, "factor")) contrasts(el))
    contr <- contr[!unlist(lapply(contr, is.null))]
    X <- model.matrix(model, X)
    if(ncol(X) == 0L) stop("no coefficients to fit")
    y <- eval(model[[2L]], dataMod)
    N <- nrow(X)
    p <- ncol(X)				# number of coefficients
    parAssign <- attr(X, "assign")
    fTerms <- terms(as.formula(model), data=data)
    namTerms <- attr(fTerms, "term.labels")
    if (attr(fTerms, "intercept") > 0) {
        namTerms <- c("(Intercept)", namTerms)
    }
    namTerms <- factor(parAssign, labels = namTerms)
    parAssign <- split(order(parAssign), namTerms)
    fixedSigma <- (controlvals$sigma > 0) ## 17-11-2015; Fixed sigma patch
    ## creating the condensed linear model
    attr(glsSt, "conLin") <-
        list(Xy = array(c(X, y), c(N, ncol(X) + 1L), list(row.names(dataMod),
	     c(colnames(X), deparse(model[[2]])))),
             dims = list(N = N, p = p, REML = as.integer(REML)), logLik = 0,
             ## 17-11-2015; Fixed sigma patch; SH Heisterkamp; Quantitative Solutions
             sigma = controlvals$sigma, fixedSigma = fixedSigma)

    ## initialization
    glsEstControl <- controlvals["singular.ok"]
    glsSt <- Initialize(glsSt, dataMod, glsEstControl)
    parMap <- attr(glsSt, "pmap")

    ##
    ## getting the fitted object, possibly iterating for variance functions
    ##
    numIter <- numIter0 <- 0L
    repeat {
        oldPars <- c(attr(glsSt, "glsFit")[["beta"]], coef(glsSt))
        if (length(coef(glsSt))) {		# needs ms()
            optRes <- if (controlvals$opt == "nlminb") {
                nlminb(c(coef(glsSt)),
                       function(glsPars) -logLik(glsSt, glsPars),
                       control = list(trace = controlvals$msVerbose,
                       iter.max = controlvals$msMaxIter))
            } else {
                optim(c(coef(glsSt)),
                      function(glsPars) -logLik(glsSt, glsPars),
                      method = controlvals$optimMethod,
                      control = list(trace = controlvals$msVerbose,
                      maxit = controlvals$msMaxIter,
                      reltol = if(numIter == 0L) controlvals$msTol
                      else 100*.Machine$double.eps))
            }
            coef(glsSt) <- optRes$par
        } else {
            optRes <- list(convergence = 0)
        }
        attr(glsSt, "glsFit") <- glsEstimate(glsSt, control = glsEstControl)
        ## checking if any updating is needed
        if (!needUpdate(glsSt)) {
            if (optRes$convergence)
                stop(optRes$message)
            break
        }
        ## updating the fit information
        numIter <- numIter + 1L
        glsSt <- update(glsSt, dataMod)
        ## calculating the convergence criterion
        aConv <- c(attr(glsSt, "glsFit")[["beta"]], coef(glsSt))
        conv <- abs((oldPars - aConv)/ifelse(aConv == 0, 1, aConv))
        aConv <- c("beta" = max(conv[1:p]))
        conv <- conv[-(1:p)]
        for(i in names(glsSt)) {
            if (any(parMap[,i])) {
                aConv <- c(aConv, max(conv[parMap[,i]]))
                names(aConv)[length(aConv)] <- i
            }
        }
        if (verbose) {
            cat("\nIteration:",numIter)
            cat("\nObjective:", format(optRes$value), "\n")
            print(glsSt)
            cat("\nConvergence:\n")
            print(aConv)
        }
        if (max(aConv) <= controlvals$tolerance) {
            break
        }
        if (numIter > controlvals$maxIter) {
            stop("maximum number of iterations reached without convergence")
        }
    }
    ## wrapping up
    glsFit <- attr(glsSt, "glsFit")
    namBeta <- names(glsFit$beta)
    ## 17-11-2015; Fixed sigma patch; SH Heisterkamp; Quantitative Solutions
    attr(glsSt, "fixedSigma") <- fixedSigma
    attr(parAssign, "varBetaFact") <- varBeta <-
        glsFit$sigma * glsFit$varBeta * sqrt((N - REML * p)/(N - p))
    varBeta <- crossprod(varBeta)
    dimnames(varBeta) <- list(namBeta, namBeta)
    ##
    ## fitted.values and residuals (in original order)
    ##
    Fitted <- fitted(glsSt)
    ## putting groups back in original order, if present
    if (!is.null(grps)) {
        grps <- grps[revOrder]
        Fitted <- Fitted[revOrder]
        Resid <- y[revOrder] - Fitted
        attr(Resid, "std") <- glsFit$sigma/varWeights(glsSt)[revOrder]
    } else {
        Resid <- y - Fitted
        attr(Resid, "std") <- glsFit$sigma/varWeights(glsSt)
    }
    names(Resid) <- names(Fitted) <- origOrder

    ## getting the approximate var-cov of the parameters
    apVar <-
        if (controlvals$apVar)
            glsApVar(glsSt, glsFit$sigma,
                     .relStep = controlvals[[".relStep"]],
                     minAbsPar = controlvals[["minAbsParApVar"]],
                     natural = controlvals[["natural"]])
        else
            "Approximate variance-covariance matrix not available"
    ## getting rid of condensed linear model and fit
    dims <- attr(glsSt, "conLin")[["dims"]]
    dims[["p"]] <- p
    attr(glsSt, "conLin") <- NULL
    attr(glsSt, "glsFit") <- NULL
    attr(glsSt, "fixedSigma") <- fixedSigma ## 17-11-2015; Fixed sigma patch; ..
    grpDta <- inherits(data, "groupedData")
    ##
    ## creating the  gls object
    ##
    structure(class = "gls",
              list(modelStruct = glsSt,
                   dims = dims,
                   contrasts = contr,
                   coefficients = glsFit[["beta"]],
                   varBeta = varBeta,
	## 17-11-2015; Fixed sigma patch; SH Heisterkamp; Quantitative Solutions
                   sigma = if(fixedSigma) controlvals$sigma else glsFit$sigma,
                   apVar = apVar,
                   logLik = glsFit$logLik,
                   numIter = if (needUpdate(glsSt)) numIter else numIter0,
                   groups = grps,
                   call = Call,
                   method = method,
                   fitted = Fitted,
                   residuals = Resid,
                   parAssign = parAssign,
		   na.action = attr(dataMod, "na.action")),
	      namBetaFull = colnames(X),
	      ## saving labels and units for plots
	      units = if(grpDta) attr(data, "units"),
	      labels= if(grpDta) attr(data, "labels"))
}

### Auxiliary functions used internally in gls and its methods

glsApVar.fullGlsLogLik <- function(Pars, object, conLin, dims, N)
{
    fixedSigma <- attr(object, "fixedSigma")
    ## logLik as a function of sigma and coef(glsSt)
    npar <- length(Pars)
    if (!fixedSigma) {
        lsigma <- Pars[npar]              # within-group std. dev.
        Pars <- Pars[-npar]
        sigma <- 0
    } else {
        sigma <- conLin$sigma
    }
    coef(object) <- Pars
    conLin <- recalc(object, conLin)
    val <- .C(gls_loglik,
              as.double(conLin$Xy),
              as.integer(unlist(dims)),
              logLik = double(1L),
              lRSS = double(1L), sigma = as.double(sigma), NAOK = TRUE)[c("logLik", "lRSS")]
    if (!fixedSigma) {
        aux <- 2 * (val[["lRSS"]] - lsigma)
        conLin[["logLik"]] + val[["logLik"]] + (N * aux - exp(aux))/2
    } else {
        val[["logLik"]]
    }
}

glsApVar <-
    function(glsSt, sigma, conLin = attr(glsSt, "conLin"),
             .relStep = .Machine$double.eps^(1/3), minAbsPar = 0,
             natural = TRUE)
{
    fixedSigma <- attr(glsSt, "fixedSigma")
    ## calculate approximate variance-covariance matrix of all parameters
    ## except the coefficients
    if (length(glsCoef <- coef(glsSt)) > 0L) {
        cSt <- glsSt[["corStruct"]]
        if (natural && !is.null(cSt) && inherits(cSt, "corSymm")) {
            cStNatPar <- coef(cSt, unconstrained = FALSE)
            class(cSt) <- c("corNatural", "corStruct")
            coef(cSt) <- log((1 + cStNatPar)/(1 - cStNatPar))
            glsSt[["corStruct"]] <- cSt
            glsCoef <- coef(glsSt)
        }
        dims <- conLin$dims
        N <- dims$N - dims$REML * dims$p
        conLin[["logLik"]] <- 0               # making sure
        Pars <- if(fixedSigma) glsCoef else c(glsCoef, lSigma = log(sigma))
        val <- fdHess(Pars, glsApVar.fullGlsLogLik, glsSt, conLin, dims, N,
                      .relStep = .relStep, minAbsPar = minAbsPar)[["Hessian"]]
        if (all(eigen(val, only.values=TRUE)$values < 0)) {
            ## negative definite - OK
            val <- solve(-val)
            nP <- names(Pars)
            dimnames(val) <- list(nP, nP)
            attr(val, "Pars") <- Pars
            attr(val, "natural") <- natural
            val
        } else {
            ## problem - solution is not a maximum
            "Non-positive definite approximate variance-covariance"
        }
    } # else NULL
}

glsEstimate <-
    function(object, conLin = attr(object, "conLin"),
             control = list(singular.ok = FALSE))
{
    dd <- conLin$dims
    p <- dd$p
    oXy <- conLin$Xy
    fixSig <- conLin$fixedSigma ## 17-11-2015; Fixed sigma patch; ..
    sigma <- conLin$sigma
    conLin <- recalc(object, conLin)	# updating for corStruct and varFunc
    val <- .C(gls_estimate,
              as.double(conLin$Xy),
              as.integer(unlist(dd)),
              beta = double(p),
              sigma = as.double(sigma), ## 17-11-2015; Fixed sigma patch; ..
              logLik = double(1L),
              varBeta = double(p * p),
              rank = integer(1),
              pivot = as.integer(1:(p + 1L)),
              NAOK = TRUE)[c("beta","sigma","logLik","varBeta", "rank", "pivot")]
    rnk <- val[["rank"]]
    rnkm1 <- rnk - 1
    if (!control$singular.ok && rnkm1 < p) {
        stop(gettextf("computed \"gls\" fit is singular, rank %s", rnk),
             domain = NA)
    }
    N <- dd$N - dd$REML * p
    namCoef <- colnames(oXy)[val[["pivot"]][1:rnkm1] + 1L]	# coef names
    varBeta <- t(array(val[["varBeta"]], c(rnkm1, rnkm1),
                       list(namCoef, namCoef)))
    beta <- val[["beta"]][1:rnkm1]
    names(beta) <- namCoef
    fitted <- c(oXy[, namCoef, drop = FALSE] %*% beta)
    resid <- oXy[, p + 1] - fitted
    ll <- conLin$logLik + val[["logLik"]]
    logLik <-
        if (!fixSig) {
            (N * (logb(N) - (1 + logb(2 * pi))))/2 + ll
            ## formula 2.21 on page 70  if sigma is estimated ML formula or 2.23 page 76 with REML
        } else {
            (-N/2) * logb(2*pi)  + ll
        }
    list(logLik = logLik, beta = beta,
         sigma = val[["sigma"]], varBeta = varBeta, fitted = fitted,
         resid = resid, auxSigma = sqrt(sum((resid)^2))/sqrt(N))
}

### Methods for standard generics

ACF.gls <-
    function(object, maxLag, resType = c("pearson", "response", "normalized"),
             form = ~1, na.action = na.fail, ...)
{
    resType <- match.arg(resType)
    res <- resid(object, type = resType)
    wchRows <- NULL
    if (is.null(grps <- getGroups(object))) {
        ## check if formula defines groups
        if (!is.null(grpsF <- getGroupsFormula(form))) {
            if (is.null(data <- getData(object))) {
                ## will try to construct
                allV <- all.vars(grpsF)
                if (length(allV) > 0L) {
                    alist <- lapply(as.list(allV), as.name)
                    names(alist) <- allV
                    alist <- c(as.list(quote(data.frame)), alist)
                    mode(alist) <- "call"
                    data <- eval(alist, sys.parent(1))
                }
            }
            grps <- model.frame(grpsF, data, na.action = na.action)
            wchRows <- !is.na(match(row.names(data), row.names(grps)))
            grps <- getGroups(grps, grpsF)
        }
    }
    if (!is.null(grps)) {
        if (!is.null(wchRows)) {
            res <- res[wchRows]
        }
        res <- split(res, grps)
    } else {
        res <- list(res)
    }
    if(missing(maxLag)) {
        maxLag <- min(c(maxL <- max(lengths(res)) - 1L,
                        as.integer(10 * log10(maxL + 1))))
    }
    val <- lapply(res,
                  function(el, maxLag) {
                      N <- maxLag + 1L
                      tt <- double(N)
                      nn <- integer(N)
                      N <- min(c(N, n <- length(el)))
                      nn[1:N] <- n + 1L - 1:N
                      ## el <- el - mean(el)
                      for(i in 1:N) {
                          el1 <- el[1:(n-i+1L)]
                          el2 <- el[i:n]
                          tt[i] <- sum(el1 * el2)
                      }
                      array(c(tt,nn), c(length(tt), 2L))
                  }, maxLag = maxLag)
    val0 <- apply(sapply(val, function(x) x[,2L]), 1, sum)
    val1 <- apply(sapply(val, function(x) x[,1L]), 1, sum)/val0
    val2 <- val1/val1[1L]
    structure(data.frame(lag = 0:maxLag, ACF = val2),
	      n.used = val0,
	      class = c("ACF", "data.frame"))
}

anova.gls <-
    function(object, ..., test = TRUE, type = c("sequential", "marginal"),
             adjustSigma = NA, Terms, L, verbose = FALSE)
{
    fixSig <- attr(object$modelStruct, "fixedSigma")
    fixSig <- !is.null(fixSig) && fixSig
    Lmiss <- missing(L)
    ## returns the likelihood ratio statistics, the AIC, and the BIC
    dots <- list(...)
    if ((rt <- length(dots) + 1L) == 1L) {    ## just one object
        if (!inherits(object,"gls"))
            stop("object must inherit from class \"gls\"")
	if(is.na(adjustSigma))
	    ## REML correction already applied to gnls objects
	    adjustSigma <- inherits(object, "gnls")
        dims <- object$dims
        N <- dims$N
        p <- dims$p
        REML <- dims$REML
        assign <- object$parAssign
        vBeta <- attr(assign, "varBetaFact")
	if (!REML && adjustSigma)
	    ## using REML-like estimate of sigma under ML
	    vBeta <- sqrt((N - p)/N) * vBeta
        c0 <- solve(t(vBeta), coef(object))
        nTerms <- length(assign)
        dDF <- N - p
        lab <- paste("Denom. DF:", dDF,"\n")
        if (missing(Terms) && Lmiss) {
            ## returns the F.table (Wald) for the fixed effects
            type <- match.arg(type)
            Fval <- Pval <- double(nTerms)
            nDF <- integer(nTerms)
            for(i in 1:nTerms) {
                nDF[i] <- length(assign[[i]])
		c0i <-
		    if (type == "sequential") # type I SS
			c0[assign[[i]]]
		    else ## "marginal"
			qr.qty(qr(vBeta[, assign[[i]], drop = FALSE]), c0)[1:nDF[i]]
                Fval[i] <- sum(c0i^2)/nDF[i]
                Pval[i] <- pf(Fval[i], nDF[i], dDF, lower.tail=FALSE)
            }
            ##
            ## fixed effects F-values, df, and p-values
            ##
	    aod <- data.frame(numDF = nDF, "F-value" = Fval, "p-value" = Pval,
			check.names=FALSE)
            rownames(aod) <- names(assign)
        } else {
            if (Lmiss) {                 # terms is given
                if (is.numeric(Terms) && all(Terms == as.integer(Terms))) {
                    if (min(Terms) < 1 || max(Terms) > nTerms) {
                        stop(gettextf("'Terms' must be between 1 and %d",
                                      nTerms), domain = NA)
                    }
                } else {
                    if (is.character(Terms)) {
                        if (any(noMatch <- is.na(match(Terms, names(assign))))) {
                            stop(sprintf(ngettext(sum(noMatch),
                                                  "term %s not matched",
                                                  "terms %s not matched"),
                                         paste(Terms[noMatch], collapse = ", ")),
                                 domain = NA)
                        }
                    } else {
                        stop("terms can only be integers or characters")
                    }
                }
                lab <-
                    paste(lab, "F-test for:",
                          paste(names(assign[Terms]),collapse=", "),"\n")
                L <- diag(p)[unlist(assign[Terms]),,drop=FALSE]
            } else {
                L <- as.matrix(L)
                if (ncol(L) == 1L) L <- t(L)     # single linear combination
                nrowL <- nrow(L)
                ncolL <- ncol(L)
                if (ncol(L) > p) {
                    stop(sprintf(ngettext(ncol(L),
                                          "'L' must have at most %d column",
                                          "'L' must have at most %d columns"),
                                 ncol(L)), domain = NA)
                }
                dmsL1 <- rownames(L)
                L0 <- array(0, c(nrowL, p), list(NULL, names(coef(object))))
                if (is.null(dmsL2 <- colnames(L))) {
                    ## assume same order as effects
                    L0[, 1:ncolL] <- L
                } else {
                    if (any(noMatch <- is.na(match(dmsL2, colnames(L0))))) {
                        stop(sprintf(ngettext(sum(noMatch),
                                              "effect %s not matched",
                                              "effects %s not matched"),
                                     paste(dmsL2[noMatch],collapse=", ")),
                             domain = NA)
                    }
                    L0[, dmsL2] <- L
                }
                L <- L0[noZeroRowL <- as.logical((L0 != 0) %*% rep(1, p)), , drop = FALSE]
                nrowL <- nrow(L)
                noZeroColL <- as.logical(c(rep(1,nrowL) %*% (L != 0)))
                rownames(L) <- if(is.null(dmsL1)) 1:nrowL else dmsL1[noZeroRowL]
                lab <- paste(lab, "F-test for linear combination(s)\n")
            }
            nDF <- sum(svd.d(L) > 0)
            c0 <- c(qr.qty(qr(vBeta %*% t(L)), c0))[1:nDF]
            Fval <- sum(c0^2)/nDF
            Pval <- pf(Fval, nDF, dDF, lower.tail=FALSE)
            aod <- data.frame(numDF = nDF, "F-value" = Fval, "p-value" = Pval,
                              check.names=FALSE)
            if (!Lmiss) {
                attr(aod, "L") <-
                    if(nrow(L) > 1) L[, noZeroColL, drop = FALSE] else L[, noZeroColL]
            }
        }
        attr(aod, "label") <- lab
        attr(aod,"rt") <- rt
        class(aod) <- c("anova.lme", "data.frame")
        aod
    }
    ##
    ## Otherwise construct the likelihood ratio and information table
    ## objects in ... may inherit from gls, lm, lmList, and lme (for now)
    ##
    else {
        Call <- match.call()
        Call[[1]] <- quote(nlme::anova.lme)
        eval.parent(Call)
    }
}

## (This is "cut'n'paste" similar to augPred.lme() in ./lme.R -- keep in sync!)
augPred.gls <-
    function(object, primary = NULL, minimum = min(primary),
             maximum = max(primary), length.out = 51L, ...)
{
    data <- eval(object$call$data)
    if (!inherits(data, "data.frame")) {
        stop(gettextf("data in %s call must evaluate to a data frame",
                      sQuote(substitute(object))), domain = NA)
    }
    if(is.null(primary)) {
        if (!inherits(data, "groupedData")) {
            stop(gettextf("%s without \"primary\" can only be used with fits of \"groupedData\" objects",
                          sys.call()[[1L]]), domain = NA)
        }
        primary <- getCovariate(data)
        pr.var <- getCovariateFormula(data)[[2L]]
    } else {
        pr.var <- asOneSidedFormula(primary)[[2L]]
        primary <- eval(pr.var, data)
    }
    prName <- c_deparse(pr.var)
    newprimary <- seq(from = minimum, to = maximum, length.out = length.out)
    groups <- getGroups(object) # much simpler here than in augPred.lme
    grName <- ".groups"
    value <- if (noGrp <- is.null(groups)) {		# no groups used
		 groups <- rep("1", length(primary))
		 data.frame(newprimary, rep("1", length(newprimary)))
	     } else {
		 ugroups <- unique(groups)
		 data.frame(rep(newprimary, length(ugroups)),
			    rep(ugroups, rep(length(newprimary), length(ugroups))))
	     }
    names(value) <- c(prName, grName)
    ## recovering other variables in data that may be needed for predictions
    ## varying variables will be replaced by their means
    summData <- gsummary(data, groups = groups)
    if (any(toAdd <- is.na(match(names(summData), names(value))))) {
        summData <- summData[, toAdd, drop = FALSE]
    }
    value[, names(summData)] <- summData[value[, 2L], ]
    pred <- predict(object, value)
    newvals <- cbind(value[, 1:2], pred)
    names(newvals)[3] <- respName <-
        deparse(resp.var <- getResponseFormula(object)[[2L]])
    orig <- data.frame(primary, groups, getResponse(object))
    names(orig) <- names(newvals)
    value <- rbind(orig, newvals)
    attributes(value[, 2]) <- attributes(groups)
    value[, ".type"] <- ordered(c(rep("original", nrow(data)),
                                  rep("predicted", nrow(newvals))),
                                levels = c("predicted", "original"))
    labs <- list(x = prName, y = respName)
    unts <- list(x = "", y = "")
    if(inherits(data, "groupedData")) {
	labs[names(attr(data, "labels"))] <- attr(data, "labels")
	unts[names(attr(data, "units"))] <- attr(data, "units")
    }
    subL <- list(Y = resp.var, X = pr.var, G = as.name(grName))
    structure(value, class = c("augPred", class(value)),
	      labels = labs,
	      units  = unts,
	      formula= if(noGrp)
			   eval (substitute(Y ~ X,     subL))
		       else eval(substitute(Y ~ X | G, subL)))
}

coef.gls <-
    function(object, allCoef = FALSE, ...)
{
    val <- object$coefficients
    if (allCoef) {
        namFull <- attr(object, "namBetaFull")
        if (length(val) != (lF <- length(namFull))) {
            aux <- rep(NA, lF)
            names(aux) <- namFull
            aux[names(val)] <- val
            val <- aux
        }
    }
    val
}

comparePred.gls <-
    function(object1, object2, primary = NULL,
             minimum = min(primary), maximum = max(primary),
             length.out = 51L, level = NULL, ...)
{
    if (length(level) > 1L) {
        stop("only one level allowed for predictions")
    }
    args <- list(object = object1,
                 primary = primary,
                 level = level,
                 length.out = length.out)
    if (!is.null(primary)) {
	  ## need to do this before forcing the evaluations
	  primary <- eval(asOneSidedFormula(primary)[[2]],
                         eval(object1$call$data))
        args[["minimum"]] <- minimum
        args[["maximum"]] <- maximum
    }
    val1 <- do.call(augPred, args)
    dm1 <- dim(val1)
    c1 <- deparse(substitute(object1))
    levels(val1[,4])[1] <- c1
    args[["object"]] <- object2
    val2 <- do.call(augPred, args)
    dm2 <- dim(val2)
    c2 <- deparse(substitute(object2))
    levels(val2[, 4L])[1] <- c2
    val2 <- val2[val2[, 4L] != "original", , drop = FALSE]
    names(val2) <- names(val1)

    if (dm1[1L] == dm2[1L]) {
        lv1 <- sort(levels(val1[, 2L]))
        lv2 <- sort(levels(val2[, 2L]))
        if ((length(lv1) != length(lv2)) || any(lv1 != lv2)) {
            stop(gettextf("%s and %s must have the same group levels", c1, c2),
                 domain = NA)
        }
        val <- rbind(val1[, -4L], val2[, -4L])
        val[, ".type"] <-
            ordered(c(as.character(val1[,4L]), as.character(val2[,4L])),
                    levels = c(c1, c2, "original"))
        attr(val, "formula") <- attr(val1, "formula")
    } else {				# one may have just "fixed"
        if (dm1[1L] > dm2[1L]) {
            mult <- dm1[1L] %/% dm2[1L]
            if ((length(levels(val2[, 2])) != 1L) ||
                (length(levels(val1[, 2])) != mult))
            {
                stop("wrong group levels")
            }
            val <-
                data.frame(c(val1[,1L], rep(val2[,1L], mult)), rep(val1[,1L], 2L),
                           c(val1[,3L], rep(val2[,3L], mult)),
                           ordered(c(as.character(val1[,4L]),
                                     rep(as.character(val2[,4L]), mult)),
                                   levels = c(c1, c2, "original")))
            attr(val, "formula") <- attr(val1, "formula")
        } else {
            mult <- dm2[1L] %/% dm1[1L]
            if ((length(levels(val1[, 2])) != 1L) ||
                (length(levels(val2[, 2])) != mult))
            {
                stop("wrong group levels")
            }
            val <-
                data.frame(c(rep(val1[,1L], mult), val2[,1L]), rep(val2[,1L], 2L),
                           c(rep(val1[,3L], mult), val2[,3L]),
                           ordered(c(rep(as.character(val1[,4L]), mult),
                                     as.character(val1[,4L])),
                                   levels = c(c1, c2, "original")))
            attr(val, "formula") <- attr(val2, "formula")
        }
    }
    class(val) <- c("comparePred", "augPred", class(val))
    attr(val, "labels") <- attr(val1, "labels")
    attr(val, "units") <- attr(val1, "units")
    val
}

fitted.gls <-
    function(object, ...)
{
    val <- napredict(object$na.action, object$fitted)
    lab <- "Fitted values"
    if (!is.null(aux <- attr(object, "units")$y)) {
        lab <- paste(lab, aux)
    }
    attr(val, "label") <- lab
    val
}


formula.gls <- function(x, ...) eval(x$call$model)

getGroups.gls <- function(object, form, level, data, sep) object$groups

getGroupsFormula.gls <- function(object, asList = FALSE, sep)
{
    if (!is.null(cSt <- object$modelStruct$corStruct))
        getGroupsFormula(cSt, asList)
    ## else NULL
}

getResponse.gls <-
    function(object, form)
{
    val <- resid(object) + fitted(object)
    if (is.null(lab <- attr(object, "labels")$y)) {
        lab <- deparse(getResponseFormula(object)[[2]])
    }
    if (!is.null(aux <- attr(object, "units")$y)) {
        lab <- paste(lab, aux)
    }
    attr(val, "label") <- lab
    val
}

intervals.gls <-
    function(object, level = 0.95, which = c("all", "var-cov", "coef"), ...)
{
    which <- match.arg(which)
    val <- list()
    dims <- object$dims
    ## 17-11-2015; Fixed sigma patch; SH Heisterkamp; Quantitative Solutions
    fixSig <- attr(object$modelStruct, "fixedSigma")
    fixSig <- !is.null(fixSig) && fixSig

    if (which != "var-cov") {		# coefficients included
        len <- -qt((1-level)/2, dims$N - dims$p) * sqrt(diag(object$varBeta))
        est <- coef(object)
        val[["coef"]] <-
            array(c(est - len, est, est + len), c(length(est), 3L),
                  list(names(est), c("lower", "est.", "upper")))
        attr(val[["coef"]], "label") <- "Coefficients:"
    }

    if (which != "coef") {		# variance-covariance included
        if (is.null(aV <- object$apVar)) {	# only sigma
            Nr <- if (inherits(object, "gnls")) {   #always REML-like sigma
                      dims$N - dims$p
                  } else {
                      dims$N - dims$REML * dims$p
                  }
            ## 17-11-2015; Fixed sigma patch; SH Heisterkamp; Quantitative Solutions
            if(!fixSig) {
                est <- object$sigma * sqrt(Nr)
                val[["sigma"]] <-
                    structure(c(lower = est/sqrt(qchisq((1+level)/2, Nr)),
                                "est."= object$sigma,
                                upper = est/sqrt(qchisq((1-level)/2, Nr))),
                              label = "Residual standard error:")
            } else {
                est <- 1
		val[["sigma"]] <- structure(setNames(rep(object$sigma, 3),
                                                     c("lower", "est.", "upper")),
                                            label = "Fixed Residual standard error:")
            }
        } else {
            if (is.character(aV)) {
                stop(gettextf("cannot get confidence intervals on var-cov components: %s",
                              aV), domain = NA)
            }
            len <- -qnorm((1-level)/2) * sqrt(diag(aV))
            est <- attr(aV, "Pars")
            nP <- length(est)
            glsSt <- object[["modelStruct"]]
            pmap <- attr(glsSt, "pmap")
            if (!all(whichKeep <- apply(pmap, 2, any))) {
                ## need to deleted components with fixed coefficients
                aux <- glsSt[whichKeep]
                class(aux) <- class(glsSt)
                attr(aux, "settings") <- attr(glsSt, "settings")
                attr(aux, "pmap") <- pmap[, whichKeep, drop = FALSE]
                glsSt <- aux
            }
            cSt <- glsSt[["corStruct"]]
            if (!is.null(cSt) && inherits(cSt, "corSymm") && attr(aV, "natural")) {
                ## converting to corNatural
                class(cSt) <- c("corNatural", "corStruct")
                glsSt[["corStruct"]] <- cSt
            }
            namG <- names(glsSt)
            auxVal <- vector("list", length(namG) + 1L)
            names(auxVal) <- c(namG, "sigma")
            aux <-
                array(c(est - len, est, est + len),
                      c(nP, 3), list(NULL, c("lower", "est.", "upper")))
	    auxVal[["sigma"]] <-
		structure(if(!fixSig) { # last param. = log(sigma)
			      s <- exp(aux[nP, ])
			      aux <- aux[-nP,, drop = FALSE]
			      s
			  } else
			      c(object$sigma, object$sigma, object$sigma),
			  label = "Residual standard error:")
            rownames(aux) <- names(coef(glsSt, FALSE))
            for(i in 1:3) {
                coef(glsSt) <- aux[,i]
                aux[,i] <- coef(glsSt, unconstrained = FALSE)
            }
            for(i in namG) {
                au.i <- aux[pmap[,i], , drop = FALSE]
                dimnames(au.i)[[1]] <- substring(dimnames(au.i)[[1]], nchar(i, "c") + 2L)
                attr(au.i, "label") <-
                    switch(i,
                           corStruct = "Correlation structure:",
                           varStruct = "Variance function:",
                           paste0(i, ":"))
                auxVal[[i]] <- au.i
            }
            val <- c(val, auxVal)
        }
    }
    structure(val, level = level, class = "intervals.gls")
}

logLik.gls <-
    function(object, REML, ...)
{
    ## 17-11-2015; Fixed sigma patch; SH Heisterkamp; Quantitative Solutions
    fixSig <- attr(object[["modelStruct"]], "fixedSigma")
    fixSig <- !is.null(fixSig) && fixSig
    p <- object$dims$p
    N <- object$dims$N
    Np <- N - p
    estM <- object$method
    if (missing(REML)) REML <- estM == "REML"
    val <- object[["logLik"]]
    if (REML && estM == "ML") { # have to correct logLik
        val <- val + (p * (log(2 * pi) + 1) + Np * log(1 - p/N) +
                      sum(log(abs(svd.d(object$varBeta))))) / 2
    }
    else if (!REML && (estM == "REML")) { # have to correct logLik
        val <- val - (p * (log(2*pi) + 1) + N * log(1 - p/N) +
                      sum(log(abs(svd.d(object$varBeta))))) / 2
    }
    structure(val,
              nall = N,
              nobs = N - REML * p,
              ## 17-11-2015; Fixed sigma patch; SH Heisterkamp; Quantitative Solutions
              df = p + length(coef(object[["modelStruct"]])) + as.integer(!fixSig),
              class = "logLik")
}

nobs.gls <- function(object, ...) object$dims$N


plot.gls <-
    function(x, form = resid(., type = "pearson") ~ fitted(.), abline,
             id = NULL, idLabels = NULL, idResType = c("pearson", "normalized"),
             grid, ...)
    ## Diagnostic plots based on residuals and/or fitted values
{
    do.call(plot.lme, as.list(match.call()[-1]))
}

predict.gls <-
    function(object, newdata, na.action = na.fail, ...)
{
    ##
    ## method for predict() designed for objects inheriting from class gls
    ##
    if (missing(newdata)) {		# will return fitted values
        return(fitted(object))
    }
    form <- getCovariateFormula(object)
    mfArgs <- list(formula = form, data = newdata, na.action = na.action)
    mfArgs$drop.unused.levels <- TRUE
    dataMod <- do.call(model.frame, mfArgs)
    ## making sure factor levels are the same as in contrasts
    contr <- object$contrasts
    for(i in names(dataMod)) {
        if (inherits(dataMod[,i], "factor") && !is.null(contr[[i]])) {
            levs <- levels(dataMod[,i])
            levsC <- dimnames(contr[[i]])[[1]]
            if (any(wch <- is.na(match(levs, levsC)))) {
                stop(sprintf(ngettext(sum(wch),
                                      "level %s not allowed for %s",
                                      "levels %s not allowed for %s"),
                             paste(levs[wch], collapse = ",")),
                     domain = NA)
            }
            attr(dataMod[,i], "contrasts") <- contr[[i]][levs, , drop = FALSE]
            ##      if (length(levs) < length(levsC)) {
            ##        if (inherits(dataMod[,i], "ordered")) {
            ##          dataMod[,i] <- ordered(as.character(dataMod[,i]), levels = levsC)
            ##        } else {
            ##          dataMod[,i] <- factor(as.character(dataMod[,i]), levels = levsC)
            ##        }
            ##      }
        }
    }
    N <- nrow(dataMod)
    if (length(all.vars(form)) > 0) {
        ##    X <- model.matrix(form, dataMod, contr)
        X <- model.matrix(form, dataMod)
    } else {
        X <- array(1, c(N, 1), list(row.names(dataMod), "(Intercept)"))
    }
    cf <- coef(object)
    val <- c(X[, names(cf), drop = FALSE] %*% cf)
    lab <- "Predicted values"
    if (!is.null(aux <- attr(object, "units")$y)) {
        lab <- paste(lab, aux)
    }
    structure(val, label = lab)
}

print.intervals.gls <-
    function(x, ...)
{
    cat(paste0("Approximate ", attr(x, "level") * 100,
               "% confidence intervals\n"))
    for(i in names(x)) {
        aux <- x[[i]]
        cat("\n ",attr(aux, "label"), "\n", sep = "")
        print(if(i == "sigma") c(aux) else as.matrix(aux), ...)
    }
    invisible(x)
}

print.gls <-
    ## method for print() used for gls objects
    function(x, ...)
{
    ## 17-11-2015; Fixed sigma patch; SH Heisterkamp; Quantitative Solutions
    fixSig <- attr(x[["modelStruct"]], "fixedSigma")
    fixSig <- !is.null(fixSig) && fixSig
    dd <- x$dims
    mCall <- x$call
    if (inherits(x, "gnls")) {
        cat("Generalized nonlinear least squares fit\n")
    } else {
        cat("Generalized least squares fit by ")
        cat(if(x$method == "REML") "REML\n" else "maximum likelihood\n")
    }
    cat("  Model:", deparse(mCall$model), "\n")
    cat("  Data:", deparse( mCall$data ), "\n")
    if (!is.null(mCall$subset)) {
        cat("  Subset:", deparse(asOneSidedFormula(mCall$subset)[[2]]),"\n")
    }
    if (inherits(x, "gnls")) {
        cat("  Log-likelihood: ", format(x$logLik), "\n", sep = "")
    } else {
        cat("  Log-", if(x$method == "REML") "restricted-" else "",
            "likelihood: ", format(x$logLik), "\n", sep = "")
    }
    cat("\nCoefficients:\n")
    print(coef(x), ...)
    cat("\n")
    if (length(x$modelStruct) > 0L) {
        print(summary(x$modelStruct), ...)
    }
    cat("Degrees of freedom:", dd[["N"]],"total;",dd[["N"]] - dd[["p"]],
        "residual\n")
    cat("Residual standard error:", format(x$sigma),"\n")
    invisible(x)
}

print.summary.gls <-
    function(x, verbose = FALSE, digits = .Options$digits, ...)
{
    dd <- x$dims
    fixSig <- attr(x[["modelStruct"]], "fixedSigma")
    fixSig <- !is.null(fixSig) && fixSig
    verbose <- verbose || attr(x, "verbose")
    mCall <- x$call
    if (inherits(x, "gnls")) {
        cat("Generalized nonlinear least squares fit\n")
    } else {
        cat("Generalized least squares fit by ")
        cat(if(x$method == "REML") "REML\n" else "maximum likelihood\n")
    }
    cat("  Model:", deparse(mCall$model), "\n")
    cat("  Data:", deparse( mCall$data ), "\n")
    if (!is.null(mCall$subset)) {
        cat("  Subset:", deparse(asOneSidedFormula(mCall$subset)[[2]]),"\n")
    }
    print(data.frame(AIC=x$AIC, BIC=x$BIC, logLik=as.vector(x$logLik),
		     row.names = " "),
	  ...)
    if (verbose) { cat("Convergence at iteration:",x$numIter,"\n") }
    if (length(x$modelStruct)) {
        cat("\n")
        print(summary(x$modelStruct), ...)
    }
    cat("\nCoefficients:\n")
    xtTab <- as.data.frame(x$tTable)
    wchPval <- match("p-value", names(xtTab))
    for(i in names(xtTab)[-wchPval]) {
        xtTab[, i] <- format(zapsmall(xtTab[, i]))
    }
    xtTab[,wchPval] <- format(round(xtTab[,wchPval], 4L))
    if (any(wchLv <- (as.double(levels(xtTab[, wchPval])) == 0))) {
        levels(xtTab[, wchPval])[wchLv] <- "<.0001"
    }
    row.names(xtTab) <- dimnames(x$tTable)[[1]]
    print(xtTab, ...)
    if (nrow(x$tTable) > 1L) {
        corr <- x$corBeta
        class(corr) <- "correlation"
        print(corr, title = "\n Correlation:", ...)
    }
    cat("\nStandardized residuals:\n")
    print(x$residuals, ...)
    cat("\n")
    cat("Residual standard error:", format(x$sigma),"\n")
    cat("Degrees of freedom:", dd[["N"]],"total;", dd[["N"]] - dd[["p"]],
        "residual\n")
    invisible(x)
}

residuals.gls <-
    function(object, type = c("response", "pearson", "normalized"), ...)
{
    type <- match.arg(type)
    val <- object$residuals
    if (type != "response") {
        val <- val/attr(val, "std")
        lab <- "Standardized residuals"
        if (type == "normalized") {
            if (!is.null(cSt <- object$modelStruct$corStruct)) {
                ## normalize according to inv-trans factor
                val <- recalc(cSt, list(Xy = as.matrix(val)))$Xy[, 1]
                lab <- "Normalized residuals"
            }
        }
    } else {
        lab <- "Residuals"
        if (!is.null(aux <- attr(object, "units")$y))
            lab <- paste(lab, aux)
    }
    attr(val, "label") <- lab
    if (!is.null(object$na.action)) {
        res <- naresid(object$na.action, val)
        attr(res, "std") <- naresid(object$na.action, attr(val, "std"))
        attr(res, "label") <- attr(val, "label")
        res
    } else val
}

summary.gls <- function(object, verbose = FALSE, ...) {
    ## 17-11-2015; Fixed sigma patch; SH Heisterkamp; Quantitative Solutions
    fixSig <- attr(object[["modelStruct"]], "fixedSigma")
    fixSig <- !is.null(fixSig) && fixSig
    ##
    ## generates an object used in the print.summary method for lme
    ##
    ##  variance-covariance estimates for coefficients
    ##
    stdBeta <- sqrt(diag(as.matrix(object$varBeta)))
    corBeta <- t(object$varBeta/stdBeta)/stdBeta
    ##
    ## coefficients, std. deviations and z-ratios
    ##
    beta <- coef(object)
    dims <- object$dims
    dimnames(corBeta) <- list(names(beta),names(beta))
    object$corBeta <- corBeta
    tTable <- data.frame(beta, stdBeta, beta/stdBeta, beta)
    dimnames(tTable)<-
        list(names(beta),c("Value","Std.Error","t-value","p-value"))
    tTable[, "p-value"] <- 2 * pt(-abs(tTable[,"t-value"]), dims$N - dims$p)
    object$tTable <- as.matrix(tTable)
    ##
    ## residuals
    ##
    resd <- resid(object, type = "pearson")
    if (length(resd) > 5) {
        resd <- quantile(resd, na.rm = TRUE)
        names(resd) <- c("Min","Q1","Med","Q3","Max")
    }
    object$residuals <- resd
    ##
    ## generating the final object
    ##
    aux <- logLik(object)
    structure(c(object, list(BIC = BIC(aux), AIC = AIC(aux))),
	      verbose = verbose,
	      class = c("summary.gls", class(object)))
}

update.gls <-
    function (object, model., ..., evaluate = TRUE)
{
    call <- object$call
    if (is.null(call))
	stop("need an object with call component")
    extras <- match.call(expand.dots = FALSE)$...
    if (!missing(model.))
	call$model <- update.formula(formula(object), model.)
    if(length(extras) > 0L) {
        ## the next two lines allow partial matching of argument names
        ## in the update.  This is nonstandard but required for examples
        ## in chapter 5 of Pinheiro and Bates (2000).
        glsa <- names(as.list(args(gls)))
        names(extras) <- glsa[pmatch(names(extras), glsa[-length(glsa)])]
	existing <- !is.na(match(names(extras), names(call)))
	## do these individually to allow NULL to remove entries.
	for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
	if(any(!existing)) {
	    call <- c(as.list(call), extras[!existing])
	    call <- as.call(call)
	}
    }
    if(evaluate) eval(call, parent.frame())
    else call
}

Variogram.gls <-
    function(object, distance, form = ~1,
             resType = c("pearson", "response", "normalized"),
             data, na.action = na.fail, maxDist, length.out = 50,
             collapse = c("quantiles", "fixed", "none"), nint = 20, breaks,
             robust = FALSE, metric = c("euclidean", "maximum", "manhattan"),
             ...)
{
    resType <- match.arg(resType)
    ## checking if object has a corSpatial element
    csT <- object$modelStruct$corStruct
    wchRows <- NULL
    if (missing(distance)) {
        if (missing(form) && inherits(csT, "corSpatial")) {
            distance <- getCovariate(csT)
            grps <- getGroups(object)
        } else {
            metric <- match.arg(metric)
            if (missing(data)) {
                data <- getData(object)
            }
            if (is.null(data)) {			# will try to construct
                allV <- all.vars(form)
                if (length(allV) > 0L) {
                    alist <- lapply(as.list(allV), as.name)
                    names(alist) <- allV
                    alist <- c(as.list(quote(data.frame)), alist)
                    mode(alist) <- "call"
                    data <- eval(alist, sys.parent(1))
                }
            }
            grpsF <- getGroupsFormula(form)
            grps <- NULL
            if (is.null(grpsF) || is.null(grps <- getGroups(data, grpsF))) {
                ## try to get from object
                grps <- getGroups(object)
            }
            covForm <- getCovariateFormula(form)
            covar <-
                if (length(all.vars(covForm)) > 0L) {
                    if (attr(terms(covForm), "intercept") == 1L) {
                        covForm <- eval(substitute(~ CV - 1, list(CV = covForm[[2]])))
                    }
                    covar <- model.frame(covForm, data, na.action = na.action)
                    ## making sure grps is consistent
                    wchRows <- !is.na(match(row.names(data), row.names(covar)))
                    if (!is.null(grps)) {
                        grps <- grps[wchRows, drop = TRUE]
                    }
                    as.data.frame(unclass(model.matrix(covForm, covar)))
                } else if (is.null(grps))
                    1:nrow(data)
                else
                    data.frame(dist = unlist(tapply(rep(1, nrow(data)), grps, cumsum)))

            distance <-
                if (is.null(grps))
                    dist(as.matrix(covar), method = metric)
                else {
                    covar <- split(covar, grps)
                    ## getting rid of 1-observation groups
                    covar <- covar[sapply(covar, function(el) nrow(as.matrix(el))) > 1]
                    lapply(covar,
                           function(el, metric) dist(as.matrix(el), method=metric),
                           metric = metric)
                }
        }
    }
    res <- resid(object, type = resType)
    if (!is.null(wchRows)) {
        res <- res[wchRows]
    }
    if (is.null(grps)) {
        val <- Variogram(res, distance)
    } else {
        res <- split(res, grps)
        res <- res[lengths(res) > 1L] # no 1-observation groups
        levGrps <- levels(grps)
        val <- structure(vector("list", length(levGrps)), names = levGrps)
        for(i in levGrps) {
            val[[i]] <- Variogram(res[[i]], distance[[i]])
        }
        val <- do.call(rbind, val)
    }
    if (!missing(maxDist)) {
        val <- val[val$dist <= maxDist, ]
    }
    collapse <- match.arg(collapse)
    if (collapse != "none") {             # will collapse values
        dst <- val$dist
        udist <- sort(unique(dst))
        ludist <- length(udist)
        if (!missing(breaks)) {
            if (min(breaks) > udist[1L]) {
                breaks <- c(udist[1L], breaks)
            }
            if (max(breaks) < udist[2L]) {
                breaks <- c(breaks, udist[2L])
            }
            if (!missing(nint) && nint != (length(breaks) - 1L)) {
                stop("'nint' is not consistent with 'breaks'")
            }
            nint <- length(breaks) - 1L
        }
        if (nint < ludist) {
            if (missing(breaks)) {
                if (collapse == "quantiles") {    # break into equal groups
                    breaks <- unique(quantile(dst, seq(0, 1, 1/nint)))
                } else {                          # fixed length intervals
                    breaks <- seq(udist[1L], udist[length(udist)], length = nint + 1)
                }
            }
            cutDist <- cut(dst, breaks)
        } else {
            cutDist <- dst
        }
        val <- lapply(split(val, cutDist),
                      function(el, robust) {
                          nh <- nrow(el)
                          vrg <- el$variog
                          if (robust) {
                              vrg <- ((mean(vrg^0.25))^4)/(0.457+0.494/nh)
                          } else {
                              vrg <- mean(vrg)
                          }
                          dst <- median(el$dist)
                          data.frame(variog = vrg, dist = dst)
                      }, robust = robust)
        val <- do.call(rbind, as.list(val))
        val$n.pairs <- unclass(table(na.omit(cutDist)))
    }
    row.names(val) <- 1:nrow(val)
    if (inherits(csT, "corSpatial") && resType != "normalized") {
        ## will keep model variogram
        sig2 <- if (resType == "pearson") 1 else object$sigma^2
        attr(val, "modelVariog") <-
            Variogram(csT, sig2 = sig2, length.out = length.out)
    }
    structure(val, collapse = collapse != "none",
	      class = c("Variogram", "data.frame"))
}

###*### glsStruct - a model structure for gls fits

glsStruct <-
    ## constructor for glsStruct objects
    function(corStruct = NULL, varStruct = NULL)
{
    val <- list(corStruct = corStruct, varStruct = varStruct)
    val <- val[!sapply(val, is.null)]	# removing NULL components
    class(val) <- c("glsStruct", "modelStruct")
    val
}

##*## glsStruct methods for standard generics

fitted.glsStruct <-
    function(object, glsFit = attr(object, "glsFit"), ...)
{
    glsFit[["fitted"]]
}

Initialize.glsStruct <-
    function(object, data, control = list(singular.ok = FALSE), ...)
{
    if (length(object)) {
        object[] <- lapply(object, Initialize, data)
        theta <- lapply(object, coef)
        len <- lengths(theta)
        num <- seq_along(len)
        pmap <-
            if (sum(len) > 0)
                outer(rep(num, len), num, "==")
            else
                array(FALSE, c(1, length(len)))
        dimnames(pmap) <- list(NULL, names(object))
        attr(object, "pmap") <- pmap
        attr(object, "glsFit") <-
            glsEstimate(object, control = control)
        if (needUpdate(object)) {
            object <- update(object, data)
        }
    }
    object
}

logLik.glsStruct <-
    function(object, Pars, conLin = attr(object, "conLin"), ...)
{
    coef(object) <- Pars			# updating parameter values
    conLin <- recalc(object, conLin)	# updating conLin
    val <- .C(gls_loglik,
              as.double(conLin[["Xy"]]),
              as.integer(unlist(conLin[["dims"]])),
              logLik = as.double(conLin[["logLik"]]),
              ## 17-11-2015; Fixed sigma patch; SH Heisterkamp; Quantitative Solutions
              double(1L), as.double(conLin$sigma), NAOK = TRUE)
    val[["logLik"]]
}

residuals.glsStruct <-
    function(object, glsFit = attr(object, "glsFit"), ...)
{
    glsFit[["resid"]]
}

varWeights.glsStruct <-
    function(object)
{
    if (is.null(object$varStruct)) rep(1, attr(object, "conLin")$dims$N)
    else varWeights(object$varStruct)
}

## Auxiliary control functions

glsControl <-
    ## Control parameters for gls
    function(maxIter = 50L, msMaxIter = 200L, tolerance = 1e-6, msTol = 1e-7,
             msVerbose = FALSE, singular.ok = FALSE,
             returnObject = FALSE,
             apVar = TRUE, .relStep = .Machine$double.eps^(1/3),
	     opt = c("nlminb", "optim"),  optimMethod = "BFGS",
             minAbsParApVar = 0.05, natural = TRUE, sigma = NULL)
{
    ## 17-11-2015; Fixed sigma patch; SH Heisterkamp; Quantitative Solutions
    if(is.null(sigma))
	sigma <- 0
    else {
	if(!is.finite(sigma) || length(sigma) != 1 || sigma < 0)
	    stop("Within-group std. dev. must be a positive numeric value")
	## if(missing(apVar)) apVar <- FALSE # not yet implemented
    }
    list(maxIter = maxIter, msMaxIter = msMaxIter, tolerance = tolerance,
         msTol = msTol, msVerbose = msVerbose,
         singular.ok = singular.ok,
         returnObject = returnObject, apVar = apVar,
         minAbsParApVar = minAbsParApVar, .relStep = .relStep,
         opt = match.arg(opt),
	 optimMethod = optimMethod, natural = natural, sigma = sigma)
}
