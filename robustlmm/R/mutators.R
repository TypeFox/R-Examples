### Set sigma
setSigma <- function(object, value) {
    object@pp$setSigma(value)
    invisible(object)
}

### Keep this for backwards compatibility
`sigma<-` <- function(object, value) {
    setSigma(object, value)
    dd <- object@devcomp$dims
    object@devcomp$cmp[ifelse(dd["REML"], "sigmaREML", "sigmaML")] <- value
    object
}

### Set fixed effects
setFixef <- function(object, value, ...) {
    object@pp$beta <- value 
    ## update mu and wtres
    object@resp$updateMu(.mu(object))
    invisible(object)
}
## keep fixef<- for backwards compatibility
`fixef<-` <- function(object, value) {
    setFixef(object, value)
    object@beta <- object@pp$beta
    object
}

### Set u
setU <- function(object, value) {
    object@pp$setU(value)
    invisible(object)
}
## keep u<- for backwards compatiblity
`u<-` <- function(object, value) {
    setU(object, value)
    object@b.s <- object@pp$b.s
    object@b.r <- object@pp$b.r
   object
}

### Set b
setB <- function(object, value) {
    object@pp$setB(value)
    invisible(object)
}
## keep b<- for backwards compatibility
`b<-` <- function(object, value) {
    setB(object, value)
    object@b.s <- object@pp$b.s
    object@b.r <- object@pp$b.r
    object
}

### Set theta
setTheta <- function(object, value, eps = 1e-7, fit.effects = TRUE,
                     update.sigma = fit.effects, ...) {
    stopifnot(length(value) == len(object, "theta"))
    if (!fit.effects && update.sigma)
        stop("fit.effects == FALSE implies update.sigma == FALSE")
    ## set thetas smaller than eps to 0
    ## (otherwise numerical stability might be critical!)
    if (any(idx <- abs(value) < eps & object@lower >= 0))
        value[idx] <- 0
    stopifnot(all(value >= object@lower))
    ## fix theta
    offset <- 0
    for (trm in object@cnms) {
        nc <- length(trm)
        ## current index
        idx.cur <- seq.int(length(value)) %in% (offset + seq.int((nc * (nc + 1))/2))
        ## all covariances
        idx.cov <- object@lower >= 0 & idx.cur
        ## if all vc are also 0 set correlation = 0
        if (all(value[idx.cov] == 0))
            value[idx.cur] <- 0
        offset <- offset + length(trm)
    }
    ## cat("Setting theta to", value, "\n")
    object@pp$setTheta(value)
    ## update sigma, fit effects?
    if (fit.effects) {
        if (update.sigma) {
            updateSigma(object, fit.effects = fit.effects)
            ##cat("Updated sigma to", object@pp$sigma, "\n")
        } else {
            ## cat("Starting values:", object@pp$beta, object@pp$b.s, "\n")
            fitEffects(c(object@pp$beta, object@pp$b.s), object)
        }
    } else {
        ## update b.s according to the new Lambda
        setB(object, object@pp$b.r)
        setFixef(object, object@pp$beta)
        ## no need to update sigma here
    }
    invisible(object)
}

## keep this for backwards compatibility
`theta<-` <- function(object, ..., value) {
    setTheta(object, value, ...)
    object@theta <- theta(object)
    updateWeights(object)
}
