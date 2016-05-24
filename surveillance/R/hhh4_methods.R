################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Standard methods for hhh4-fits
###
### Copyright (C) 2010-2012 Michaela Paul, 2012-2016 Sebastian Meyer
### $Revision: 1618 $
### $Date: 2016-03-11 09:53:42 +0100 (Fre, 11. Mär 2016) $
################################################################################

## NOTE: we also apply print.hhh4 in print.summary.hhh4()
print.hhh4 <- function (x, digits = max(3, getOption("digits")-3), ...)
{
    if (!x$convergence) {
        cat('Results are not reliable! Try different starting values.\n')
        return(invisible(x))
    }
    if (!is.null(x$call)) {
        cat("\nCall: \n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
            "\n\n", sep = "")
    }
    if (x$dim["random"] > 0) {
        cat('Random effects:\n')
        .printREmat(if (is.null(x$REmat)) .getREmat(x) else x$REmat,
                    digits = digits)
        cat("\nFixed effects:\n")
    } else if (x$dim["fixed"] > 0) {
        cat("Coefficients:\n")
    }
    if (x$dim["fixed"] > 0) {
        print.default(
            format(if (is.null(x$fixef)) fixef.hhh4(x, ...) else x$fixef,
                   digits=digits),
            quote = FALSE, print.gap = 2)
    } else cat("No coefficients\n")
    cat("\n")
    invisible(x)
}

## get estimated covariance matrix of random effects
.getREmat <- function (object)
{
    ## return NULL if model has no random effects
    if (is.null(REmat <- object$Sigma)) return(NULL)

    ## hhh4()$Sigma is named since r791 only -> derive names from Sigma.orig
    if (is.null(dimnames(REmat)))
        dimnames(REmat) <- rep.int(
            list(sub("^sd\\.", "",
                     names(object$Sigma.orig)[seq_len(nrow(REmat))])), 2L)
    
    attr(REmat, "correlation") <- cov2cor(REmat)
    attr(REmat, "sd") <- sqrt(diag(REmat))
    REmat
}
.printREmat <- function (REmat, digits = 4)
{
    V <- round(diag(REmat), digits=digits)
    corr <- round(attr(REmat, "correlation"), digits=digits)
    corr[upper.tri(corr,diag=TRUE)] <- ""
    V.corr <- cbind(V, corr, deparse.level=0)
    colnames(V.corr) <- c("Var", "Corr", rep.int("", ncol(corr)-1L))
    print.default(V.corr, quote=FALSE)
}

summary.hhh4 <- function (object, maxEV = FALSE, ...)
{
    ## do not summarize results in case of non-convergence
    if (!object$convergence) {
        cat('Results are not reliable! Try different starting values.\n')
	return(invisible(object))
    }
    ret <- c(object[c("call", "convergence", "dim", "loglikelihood", "margll",
                      "lags", "nTime", "nUnit")],
             list(fixef = fixef.hhh4(object, se=TRUE, ...),
                  ranef = ranef.hhh4(object, ...),
                  REmat = .getREmat(object),
                  AIC   = AIC(object),
                  BIC   = BIC(object),
                  maxEV_range = if (maxEV) unique(range(getMaxEV(object)))))
    class(ret) <- "summary.hhh4"
    return(ret)
}

print.summary.hhh4 <- function (x, digits = max(3, getOption("digits")-3), ...)
{
    ## x$convergence is always TRUE if we have a summary
    print.hhh4(x) # also works for summary.hhh4-objects

    if (!is.null(x$maxEV_range))
        cat("Epidemic dominant eigenvalue: ",
            paste(sprintf("%.2f", x$maxEV_range), collapse = " -- "), "\n\n")
    if(x$dim["random"]==0){
        cat('Log-likelihood:  ',round(x$loglikelihood,digits=digits-2),'\n')  
        cat('AIC:             ',round(x$AIC,digits=digits-2),'\n')
        cat('BIC:             ',round(x$BIC,digits=digits-2),'\n\n')
    } else {
        cat('Penalized log-likelihood: ',round(x$loglikelihood,digits=digits-2),'\n')  
        cat('Marginal log-likelihood:  ',round(x$margll,digits=digits-2),'\n\n')        
    }
    cat('Number of units:       ', x$nUnit, '\n')
    cat('Number of time points: ', x$nTime, '\n')
    if (!is.null(x$lags)) { # only available since surveillance 1.8-0
        if (!is.na(x$lags["ar"]) && x$lags["ar"] != 1)
            cat("Non-default autoregressive lag:  ", x$lags[["ar"]], "\n")
        if (!is.na(x$lags["ne"]) && x$lags["ne"] != 1)
            cat("Non-default neighbor-driven lag: ", x$lags[["ne"]], "\n")
    }
    cat("\n")
    invisible(x)
}

terms.hhh4 <- function (x, ...)
{
    if (is.null(x$terms))
        interpretControl(x$control,x$stsObj) else x$terms
}

nobs.hhh4 <- function (object, ...) {
    if (object$convergence) object$nObs else NA_real_
}

logLik.hhh4 <- function(object, ...)
{
    val <- if (object$convergence) object$loglikelihood else {
        warning("algorithm did not converge")
        NA_real_
    }
    attr(val, "df") <- if (object$dim["random"])
        NA_integer_ else object$dim[["fixed"]]  # use "[[" to drop the name
    attr(val, "nobs") <- nobs.hhh4(object)
    class(val) <- "logLik"
    val
}

coef.hhh4 <- function(object, se=FALSE,
                      reparamPsi=TRUE, idx2Exp=NULL, amplitudeShift=FALSE, ...)
{
    if (identical(object$control$family, "Poisson")) reparamPsi <- FALSE
    coefs <- object$coefficients
    coefnames <- names(coefs)
    idx <- getCoefIdxRenamed(coefnames, reparamPsi, idx2Exp, amplitudeShift,
                             warn=!se)
    
    ## transform and rename
    if (length(idx$Psi)) {
        coefs[idx$Psi] <- exp(-coefs[idx$Psi])  # -log(overdisp) -> overdisp
        coefnames[idx$Psi] <- names(idx$Psi)
    }
    if (length(idx$toExp)) {
        coefs[idx$toExp] <- exp(coefs[idx$toExp])
        coefnames[idx$toExp] <- names(idx$toExp)
    }
    if (length(idx$AS)) {
        coefs[idx$AS] <- sinCos2amplitudeShift(coefs[idx$AS])
        coefnames[idx$AS] <- names(idx$AS)
    }
    ## set new names
    names(coefs) <- coefnames
    
    if (se) {
        cov <- vcov.hhh4(object, reparamPsi=reparamPsi, idx2Exp=idx2Exp,
                         amplitudeShift=amplitudeShift)
        cbind("Estimate"=coefs, "Std. Error"=sqrt(diag(cov)))
    } else coefs
}

vcov.hhh4 <- function (object,
                       reparamPsi=TRUE, idx2Exp=NULL, amplitudeShift=FALSE, ...)
{
    if (identical(object$control$family, "Poisson")) reparamPsi <- FALSE
    idx <- getCoefIdxRenamed(names(object$coefficients),
                             reparamPsi, idx2Exp, amplitudeShift, warn=FALSE)
    newcoefs <- coef.hhh4(object, se=FALSE, reparamPsi=reparamPsi,
                          idx2Exp=idx2Exp, amplitudeShift=amplitudeShift)

    ## Use multivariate Delta rule => D %*% vcov %*% t(D), D: Jacobian.
    ## For idx2Exp and reparamPsi, we only transform coefficients independently,
    ## i.e. D is diagonal (with elements 'd')
    d <- rep.int(1, length(newcoefs))
    if (length(idx$Psi)) # h = exp(-psi), h' = -exp(-psi)
        d[idx$Psi] <- -newcoefs[idx$Psi]
    if (length(idx$toExp)) # h = exp(coef), h' = exp(coef)
        d[idx$toExp] <- newcoefs[idx$toExp]
    ## For the amplitude/shift-transformation, D is non-diagonal
    vcov <- if (length(idx$AS)) {
        D <- diag(d, length(d))
        D[idx$AS,idx$AS] <- jacobianAmplitudeShift(newcoefs[idx$AS])
        D %*% object$cov %*% t(D)
    } else t(t(object$cov*d)*d)  # 30 times faster than via matrix products
        
    dimnames(vcov) <- list(names(newcoefs), names(newcoefs))
    vcov
}

getCoefIdxRenamed <- function (coefnames, reparamPsi=TRUE, idx2Exp=NULL,
                               amplitudeShift=FALSE, warn=TRUE)
{
    ## indexes of overdispersion parameters
    idxPsi <- if (reparamPsi) {
        idxPsi <- grep("-log(overdisp", coefnames, fixed=TRUE)
        ## change labels from "-log(overdisp.xxx)" to "overdisp.xxx"
        names(idxPsi) <- substr(coefnames[idxPsi], start=6,
                                stop=nchar(coefnames[idxPsi])-1L)
        if (length(idxPsi) == 0L) { # backward compatibility (internal psi coef
                                    # was named "overdisp" prior to r406)
            idxPsi <- grep("^overdisp", coefnames)
            names(idxPsi) <- coefnames[idxPsi]
        }
        idxPsi
    } else NULL

    ## indexes of sine-cosine coefficients
    idxAS <- if (amplitudeShift) {
        idxAS <- sort(c(grep(".sin(", coefnames, fixed=TRUE),
                        grep(".cos(", coefnames, fixed=TRUE)))
        names(idxAS) <- sub(".sin", ".A", coefnames[idxAS], fixed=TRUE)
        names(idxAS) <- sub(".cos", ".s", names(idxAS), fixed=TRUE)
        idxAS
    } else NULL

    ## indexes of coefficients to exp()-transform
    if (isTRUE(idx2Exp)) {
        idxLogCovar <- grep(".log(", coefnames, fixed = TRUE)
        idx2Exp <- setdiff(seq_along(coefnames), c(idxLogCovar, idxPsi, idxAS))
    } else if (length(idx2Exp)) {
        stopifnot(is.vector(idx2Exp, mode = "numeric"))
        ## index sets must be disjoint
        if (length(idxOverlap <- intersect(c(idxPsi, idxAS), idx2Exp))) {
            if (warn)
                warning("following 'idx2Exp' were ignored due to overlap: ",
                        paste(idxOverlap, collapse=", "))
            idx2Exp <- setdiff(idx2Exp, idxOverlap)
        }
    }
    if (length(idx2Exp))
        names(idx2Exp) <- paste0("exp(", coefnames[idx2Exp], ")")

    ## done
    list(Psi=idxPsi, AS=idxAS, toExp=idx2Exp)
}

fixef.hhh4 <- function (object,...)
{
    if (object$dim[1L] > 0) {
        head(coef.hhh4(object, ...), object$dim[1L])
    } else NULL
}

ranef.hhh4 <- function (object, tomatrix = FALSE, ...)
{
    if (object$dim[2L] > 0){
        ranefvec <- tail(coef.hhh4(object, ...), object$dim[2L])
    } else return(NULL)
    if (!tomatrix) return(ranefvec)

    ## transform to a nUnits x c matrix (c %in% 1:3)
    model <- terms.hhh4(object)
    idxRE <- model$indexRE
    idxs <- unique(idxRE)
    names(idxs) <- model$namesFE[idxs]
    mat <- sapply(idxs, function (idx) {
        RE <- ranefvec[idxRE==idx]
        Z <- model$terms["Z.intercept",][[idx]]
        "%m%" <- get(model$terms["mult",][[idx]])
        Z %m% RE
    })
    rownames(mat) <- colnames(model$response)
    return(mat)
}

## adaption of stats::confint.default authored by the R Core Team
confint.hhh4 <- function (object, parm, level = 0.95,
                          reparamPsi=TRUE, idx2Exp=NULL, amplitudeShift=FALSE,
                          ...)
{
    cf <- coef.hhh4(object, se=TRUE, reparamPsi=reparamPsi, idx2Exp=idx2Exp,
                    amplitudeShift=amplitudeShift, ...)
    ## CAVE: random intercepts have no names (all "")
    if (missing(parm))
        parm <- seq_len(nrow(cf))
    pnames <- if (is.numeric(parm)) rownames(cf)[parm] else parm
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    pct <- paste(format(100*a, trim=TRUE, scientific=FALSE, digits=3), "%")
    fac <- qnorm(a)
    ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(pnames, pct))
    ses <- cf[parm,2]
    ci[] <- cf[parm,1] + ses %o% fac
    ci
}


## mean predictions for a subset of 1:nrow(object$stsObj)
predict.hhh4 <- function(object, newSubset = object$control$subset,
                         type = "response", ...)
{
    if (type == "response" &&
        all((m <- match(newSubset, object$control$subset, nomatch=0L)) > 0)) {
        ## we can extract fitted means from object
        object$fitted.values[m,,drop=FALSE]
    } else { ## means for time points not fitted (not part of object$control$subset)
        predicted <- meanHHH(coef.hhh4(object, reparamPsi=FALSE),
                             terms.hhh4(object),
                             subset=newSubset)
        if (type=="response") predicted$mean else {
            type <- match.arg(type, names(predicted))
            predicted[[type]]
        }
    }
}


### refit hhh4-model
## ...: arguments modifying the original control list
## S: a named list to adjust the number of harmonics of the three components
## subset.upper: refit on a subset of the data up to that time point
## use.estimates: use fitted parameters as new start values

update.hhh4 <- function (object, ..., S = NULL, subset.upper = NULL,
                         use.estimates = TRUE, evaluate = TRUE)
{
    control <- object$control

    ## first modify the control list according to the components in ...
    extras <- list(...)
    control <- modifyList(control, extras)

    ## adjust start values
    control$start <- if (use.estimates) { # use parameter estimates
        hhh4coef2start(object)
    } else local({ # re-use previous 'start' specification
        ## for pre-1.8-2 "hhh4" objects,
        ## object$control$start is not necessarily a complete list:
        template <- eval(formals(hhh4)$control$start)
        template[] <- object$control$start[names(template)]
        template
    })
    ## and update according to an extra 'start' argument
    if (!is.null(extras[["start"]])) {
        if (!is.list(extras$start) || is.null(names(extras$start))) {
            stop("'start' must be a named list, see 'help(\"hhh4\")'")
        }
        control$start[] <- mapply(
            FUN = function (now, extra) {
                if (is.null(names(extra))) {
                    extra
                } else { # can retain non-extra values
                    now[names(extra)] <- extra
                    now
                }
            },
            control$start, extras$start[names(control$start)],
            SIMPLIFY = FALSE, USE.NAMES = FALSE
        )
    }
    ## update initial values of parametric weight function
    if (use.estimates && length(coefW <- coefW(object)) &&
        ! "weights" %in% names(extras$ne)) { # only if function is unchanged
        control$ne$weights$initial <- coefW
    }

    ## adjust seasonality
    if (!is.null(S)) {
        stopifnot(is.list(S), !is.null(names(S)),
                  names(S) %in% c("ar", "ne", "end"))
        control[names(S)] <- mapply(function (comp, S) {
            comp$f <- addSeason2formula(removeSeasonFromFormula(comp$f),
                                        period = object$stsObj@freq, S = S)
            comp
        }, control[names(S)], S, SIMPLIFY=FALSE, USE.NAMES=FALSE)
    }

    ## restrict fit to those epochs of control$subset which are <=subset.upper
    if (isScalar(subset.upper))
        control$subset <- control$subset[control$subset <= subset.upper]


    ## fit the updated model or just return the modified control list
    if (evaluate) { 
        hhh4(stsObj = object$stsObj, control = control)
    } else {
        control
    }
}

## remove sine-cosine terms from a formula
## f: usually a model "formula", but can generally be of any class for which
##    terms() and formula() apply
removeSeasonFromFormula <- function (f)
{
    fterms <- terms(f, keep.order = TRUE)
    ## search sine-cosine terms of the forms "sin(..." and "fe(sin(..."
    idxSinCos <- grep("^(fe\\()?(sin|cos)\\(", attr(fterms, "term.labels"))
    formula(if (length(idxSinCos)) fterms[-idxSinCos] else f)
}

## remove all temporal terms from a formula
removeTimeFromFormula <- function (f, timevar = "t") {
    fterms <- terms(f, keep.order = TRUE)
    containsTime <- vapply(attr(fterms, "variables")[-1L],
                           FUN = function (x) timevar %in% all.vars(x),
                           FUN.VALUE = TRUE, USE.NAMES = FALSE)
    formula(if (any(containsTime)) fterms[!containsTime] else f)
}


## convert fitted parameters to a list suitable for control$start
hhh4coef2start <- function (fit)
{
    res <- list(fixed = fit$coefficients[seq_len(fit$dim[1L])],
                random = fit$coefficients[fit$dim[1L]+seq_len(fit$dim[2L])],
                sd.corr = fit$Sigma.orig)
    if (any(!nzchar(names(res$random)))) { # no names pre 1.8-2
        names(res$random) <- NULL
    }
    res
}

## extract coefficients in a list
coeflist.hhh4 <- function (x, ...)
{
    ## determine number of parameters by parameter group
    model <- terms.hhh4(x)
    dim.fe.group <- unlist(model$terms["dim.fe",], recursive = FALSE, use.names = FALSE)
    dim.re.group <- unlist(model$terms["dim.re",], recursive = FALSE, use.names = FALSE)
    nFERE <- lapply(X = list(fe = dim.fe.group, re = dim.re.group),
           FUN = function (dims) {
               nParByComp <- tapply(
                   X = dims,
                   INDEX = factor(
                       unlist(model$terms["offsetComp",],
                              recursive = FALSE, use.names = FALSE),
                       levels = 1:3, labels = c("ar", "ne", "end")),
                   FUN = sum, simplify = TRUE)
               nParByComp[is.na(nParByComp)] <- 0 # component not in model 
               nParByComp
           })

    ## extract coefficients in a list (by parameter group)
    coefs <- coef.hhh4(x, se = FALSE, ...)
    list(fixed = coeflist.default(coefs[seq_len(x$dim[1L])],
             c(nFERE$fe, "neweights" = model$nd, "overdisp" = model$nOverdisp)),
         random = coeflist.default(coefs[x$dim[1L] + seq_len(x$dim[2L])],
             nFERE$re),
         sd.corr = x$Sigma.orig)
}

## extract estimated overdispersion in dnbinom() parametrization (and as matrix)
psi2size.hhh4 <- function (object, subset = object$control$subset, units = NULL)
{
    size <- sizeHHH(object$coefficients, terms.hhh4(object), subset = subset)
    if (!is.null(size) && !is.null(units)) {
        if (is.null(subset)) {
            warning("ignoring 'units' (not compatible with 'subset = NULL')")
            size
        } else {
            size[, units, drop = FALSE]
        }
    } else {
        size
    }
}

## character vector of model components that are "inModel"
componentsHHH4 <- function (object)
    names(which(sapply(object$control[c("ar", "ne", "end")], "[[", "inModel")))

## deviance residuals
residuals.hhh4 <- function (object, type = c("deviance", "response"), ...)
{
    type <- match.arg(type)
    obs <- observed(object$stsObj)[object$control$subset,]
    fit <- fitted(object)
    if (type == "response")
        return(obs - fit)

    ## deviance residuals
    ## Cf. residuals.ah, it calculates:
    ## deviance = sign(y - mean) * sqrt(2 * (distr(y) - distr(mean)))
    ## pearson = (y - mean)/sqrt(variance)
    dev.resids <- if (identical(object$control$family, "Poisson")) {
        poisson()$dev.resids
    } else {
        size <- if (identical(object$control$family, "NegBin1")) {
            psi2size.hhh4(object, subset = NULL)
        } else {
            psi2size.hhh4(object) # CAVE: a matrix -> non-standard "size"
        }
        negative.binomial(size)$dev.resids
    }

    di2 <- dev.resids(y=obs, mu=fit, wt=1)
    sign(obs-fit) * sqrt(pmax.int(di2, 0))
}

## extract the formulae of the three log-linear predictors
formula.hhh4 <- function (x, ...)
{
    lapply(x$control[c("ar", "ne", "end")], "[[", "f")
}


## decompose the fitted mean of a "hhh4" model returning an array
## with dimensions (t, i, j), where the first j index is "endemic"
decompose.hhh4 <- function (x, coefs = x$coefficients, ...)
{
    ## get three major components from meanHHH() function
    meancomps <- meanHHH(coefs, terms.hhh4(x))

    ## this contains c("endemic", "epi.own", "epi.neighbours")
    ## but we really want the mean by neighbour
    neArray <- c(meancomps$ne.exppred) * neOffsetArray(x, coefW(coefs))
    ##<- ne.exppred is (t, i) and recycled for (t, i, j)
    stopifnot(all.equal(rowSums(neArray, dims = 2), meancomps$epi.neighbours,
                        check.attributes = FALSE))
    
    ## add autoregressive part to neArray
    diagidx <- cbind(c(row(meancomps$epi.own)),
                     c(col(meancomps$epi.own)),
                     c(col(meancomps$epi.own)))
    ## usually: neArray[diagidx] == 0
    neArray[diagidx] <- neArray[diagidx] + meancomps$epi.own
    
    ## add endemic component to the array
    res <- array(c(meancomps$endemic, neArray),
          dim = dim(neArray) + c(0, 0, 1),
          dimnames = with(dimnames(neArray), list(t=t, i=i, j=c("endemic",j))))
    stopifnot(all.equal(rowSums(res, dims = 2), meancomps$mean,
                        check.attributes = FALSE))
    res
}

## get the w_{ji} Y_{j,t-1} values from a hhh4() fit
## (i.e., before summing the neighbourhood component over j)
## in an array with dimensions (t, i, j)
neOffsetArray <- function (object, pars = coefW(object),
                           subset = object$control$subset)
{
    ## initialize array ordered as (j, t, i) for apply() below
    res <- array(data = 0,
                 dim = c(object$nUnit, length(subset), object$nUnit),
                 dimnames = list(
                     "j" = colnames(object$stsObj),
                     "t" = rownames(object$stsObj)[subset],
                     "i" = colnames(object$stsObj)))
    
    ## calculate array values if the fit has an NE component
    if ("ne" %in% componentsHHH4(object)) {
        W <- getNEweights(object, pars = pars)
        Y <- observed(object$stsObj)
        tm1 <- subset - object$control$ne$lag
        is.na(tm1) <- tm1 <= 0
        tYtm1 <- t(Y[tm1,,drop=FALSE])
        res[] <- apply(W, 2L, function (wi) tYtm1 * wi)
        offset <- object$control$ne$offset
        res <- if (length(offset) > 1L) {
            offset <- offset[subset,,drop=FALSE]
            res * rep(offset, each = object$nUnit)
        } else {
            res * offset
        }
        ## stopifnot(all.equal(
        ##     colSums(res),  # sum over j
        ##     terms.hhh4(object)$offset$ne(pars)[subset,,drop=FALSE],
        ##     check.attributes = FALSE))
    }
    
    ## permute dimensions as (t, i, j)
    aperm(res, perm = c(2L, 3L, 1L), resize = TRUE)
}


## compare two hhh4 fits ignoring at least the "runtime" and "call" elements
all.equal.hhh4 <- function (target, current, ..., ignore = NULL)
{
    ignore <- unique.default(c(ignore, "runtime", "call"))
    target[ignore] <- current[ignore] <- list(NULL)
    NextMethod("all.equal")
}
