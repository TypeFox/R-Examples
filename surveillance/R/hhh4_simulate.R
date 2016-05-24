################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Simulate from a HHH4 model
###
### Copyright (C) 2012 Michaela Paul, 2013-2016 Sebastian Meyer
### $Revision: 1635 $
### $Date: 2016-03-16 12:11:48 +0100 (Mit, 16. Mär 2016) $
################################################################################


### Simulate-method for hhh4-objects

simulate.hhh4 <- function (object, # result from a call to hhh4
                           nsim=1, # number of replicates to simulate
                           seed=NULL,
                           y.start=NULL, # initial counts for epidemic components
                           subset=1:nrow(object$stsObj),
                           coefs=coef(object), # coefficients used for simulation
                           components=c("ar","ne","end"), # which comp to include
                           simplify=nsim>1, # counts array only (no full sts)
                           ...)
{
    ## Determine seed (this part is copied from stats:::simulate.lm with
    ## Copyright (C) 1995-2012 The R Core Team)
    if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)                     # initialize the RNG if necessary
    if(is.null(seed))
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    ## END seed
    
    cl <- match.call()
    theta <- if (missing(coefs)) coefs else checkCoefs(object, coefs)
    
    ## lags
    lag.ar <- object$control$ar$lag
    lag.ne <- object$control$ne$lag
    maxlag <- max(lag.ar, lag.ne)
    
    ## initial counts
    nUnits <- object$nUnit
    if (is.null(y.start)) { # set starting value to mean observed (in subset!)
        y.means <- ceiling(colMeans(observed(object$stsObj)[subset,,drop=FALSE]))
        y.start <- matrix(y.means, maxlag, nUnits, byrow=TRUE)
    } else {
        if (is.vector(y.start)) y.start <- t(y.start)
        if (ncol(y.start) != nUnits)
            stop(sQuote("y.start"), " must have nUnits=", nUnits, " columns")
        if (nrow(y.start) < maxlag)
            stop("need 'y.start' values for lag=", maxlag, " initial time points")
    }

    ## get fitted components nu_it (with offset), phi_it, lambda_it, t in subset
    model <- terms.hhh4(object)
    means <- meanHHH(theta, model, subset=subset)
    psi <- splitParams(theta,model)$overdisp

    ## weight matrix/array of the ne component
    neweights <- getNEweights(object, coefW(theta))

    ## set predictor to zero if not included ('components' argument)
    stopifnot(length(components) > 0, components %in% c("ar", "ne", "end"))
    getComp <- function (comp) {
        sel <- if (comp == "end") "endemic" else paste(comp, "exppred", sep=".")
        res <- means[[sel]]
        if (!comp %in% components) res[] <- 0
        res
    }
    ar <- getComp("ar")
    ne <- getComp("ne")
    end <- getComp("end")

    ## simulate
    simcall <- quote(
        simHHH4(ar, ne, end, psi, neweights, y.start, lag.ar, lag.ne)
        )
    if (!simplify) {
        ## result template
        res0 <- object$stsObj[subset,]
        setObserved <- function (observed) {
            res0@observed[] <- observed
            res0
        }
        simcall <- call("setObserved", simcall)
    }
    res <- if (nsim==1) eval(simcall) else
           replicate(nsim, eval(simcall),
                     simplify=if (simplify) "array" else FALSE)
    if (simplify) {
        dimnames(res)[1:2] <- list(subset, colnames(model$response))
        attr(res, "initial") <- y.start
        attr(res, "stsObserved") <- object$stsObj[subset,]
        class(res) <- "hhh4sims"
    }
    
    ## Done
    attr(res, "call") <- cl
    attr(res, "seed") <- RNGstate
    res
}


### Internal auxiliary function, which performs the actual simulation

simHHH4 <- function(ar,     # lambda_it (nTime x nUnits matrix)
                    ne,     # phi_it (nTime x nUnits matrix)
                    end,    # nu_it (nTime x nUnits matrix, offset included)
                    psi,    # overdisp param(s) or numeric(0) (psi->0 = Poisson)
                    neW,    # weight matrix/array for neighbourhood component
                    start,  # starting counts (vector of length nUnits, or
                            # matrix with nUnits columns if lag > 1)
                    lag.ar = 1,
                    lag.ne = lag.ar
                    )
{
    nTime <- nrow(end)
    nUnits <- ncol(end)

    ## simulate from Poisson or NegBin model
    rdistr <- if (length(psi)==0 ||
                  isTRUE(all.equal(psi, 0, check.attributes=FALSE))) {
        rpois
    } else {
        psi.inv <- 1/psi   # since R uses different parametrization
        ## draw 'n' samples from NegBin with mean vector 'mean' (length=nUnits)
        ## and overdispersion psi such that Variance = mean + psi*mean^2
        ## where 'size'=1/psi and length(psi) == 1 or length(mean)
        function(n, mean) rnbinom(n, mu = mean, size = psi.inv)
    }

    ## if only endemic component -> simulate independently
    if (all(ar + ne == 0)) {
        return(matrix(rdistr(length(end), end), nTime, nUnits))
    }

    ## weighted sum of counts of other (neighbouring) regions
    ## params: y - vector with (lagged) counts of regions
    ##         W - nUnits x nUnits adjacency/weight matrix (0=no neighbour)
    wSumNE <- if (is.null(neW) || all(neW == 0)) { # includes the case nUnits==1
        function (y, W) numeric(nUnits)
    } else function (y, W) .colSums(W * y, nUnits, nUnits)

    ## initialize matrices for means mu_i,t and simulated data y_i,t
    mu <- y <- matrix(0, nTime, nUnits)
    y <- rbind(start, y)
    nStart <- nrow(y) - nrow(mu)        # usually just 1 for lag=1

    ## simulate
    timeDependentWeights <- length(dim(neW)) == 3
    if (!timeDependentWeights) neWt <- neW
    for(t in seq_len(nTime)){
        if (timeDependentWeights) neWt <- neW[,,t]
        ## mean mu_i,t = lambda*y_i,t-1 + phi*sum_j wji*y_j,t-1 + nu_i,t
        mu[t,] <-
            ar[t,] * y[nStart+t-lag.ar,] +
                ne[t,] * wSumNE(y[nStart+t-lag.ne,], neWt) +
                    end[t,]
        ## Sample from Poisson/NegBin with that mean
        y[nStart+t,] <- rdistr(nUnits, mu[t,])
    }

    ## return simulated data without initial counts
    y[-seq_len(nStart),,drop=FALSE]
}


### check compatibility of a user-specified coefficient vector with model

checkCoefs <- function (object, coefs, reparamPsi=TRUE)
{
    theta <- coef(object, reparamPsi=reparamPsi)  #-> computes 1/exp(logpsi)
    if (length(coefs) != length(theta))
        stop(sQuote("coefs"), " must be of length ", length(theta))
    names(coefs) <- names(theta)
    coefs
}


### aggregate predictions over time and/or (groups of) units

aggregate.hhh4sims <- function (x, units = TRUE, time = FALSE, ..., drop = FALSE)
{
    ax <- attributes(x)

    if (time) {
        ## sum counts over the whole simulation period
        res <- colSums(x)
        ## -> a nUnits x nsim matrix -> will no longer be "hhh4sims"
        if (isTRUE(units)) { # sum over all units
            res <- colSums(res) # now a vector of length nsim
        } else if (!identical(FALSE, units)) { # sum over groups of units
            stopifnot(length(units) == dim(x)[2])
            res <- t(rowSumsBy.matrix(t(res), units))
        }
    } else {
        if (isTRUE(units)) { # sum over all units
            res <- apply(X = x, MARGIN = c(1L, 3L), FUN = sum)
            if (!drop) {
                ## restore unit dimension conforming to "hhh4sims" class
                dim(res) <- c(ax$dim[1L], 1L, ax$dim[3L])
                dnres <- ax$dimnames
                dnres[2L] <- list(NULL)
                dimnames(res) <- dnres
                ## restore attributes
                attr(res, "initial") <- as.matrix(rowSums(ax$initial))
                attr(res, "stsObserved") <- aggregate(ax$stsObserved, by = "unit")
                class(res) <- "hhh4sims"
            }
        } else if (!identical(FALSE, units)) { # sum over groups of units
            stopifnot(length(units) == dim(x)[2])
            groupnames <- names(split.default(seq_along(units), units))
            res <- apply(X = x, MARGIN = 3L, FUN = rowSumsBy.matrix, by = units)
            dim(res) <- c(ax$dim[1L], length(groupnames), ax$dim[3L])
            dnres <- ax$dimnames
            dnres[2L] <- list(groupnames)
            dimnames(res) <- dnres
            if (!drop) {
                ## restore attributes
                attr(res, "initial") <- rowSumsBy.matrix(ax$initial, units)
                attr(res, "stsObserved") <- rowSumsBy.sts(ax$stsObserved, units)
                class(res) <- "hhh4sims"
            }
        } else {
            return(x)
        }
    }
    
    ## done
    res
}

rowSumsBy.matrix <- function (x, by, na.rm = FALSE)
{
    dn <- dim(x)
    res <- vapply(X = split.default(x = seq_len(dn[2L]), f = by),
                  FUN = function (idxg)
                      .rowSums(x[, idxg, drop = FALSE],
                               dn[1L], length(idxg), na.rm = na.rm),
                  FUN.VALUE = numeric(dn[1L]), USE.NAMES = TRUE)
    if (dn[1L] == 1L) t(res) else res
}

rowSumsBy.sts <- function (x, by, na.rm = FALSE)
{
    ## map, neighbourhood, upperbound, control get lost by aggregation of units
    .sts(epoch = x@epoch, freq = x@freq, start = x@start,
        observed = rowSumsBy.matrix(x@observed, by, na.rm),
        state = rowSumsBy.matrix(x@state, by, na.rm) > 0,
        alarm = rowSumsBy.matrix(x@alarm, by, na.rm) > 0,
        populationFrac = rowSumsBy.matrix(x@populationFrac, by, na.rm),
        epochAsDate = x@epochAsDate, multinomialTS = x@multinomialTS)
}
