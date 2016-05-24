############################################################################################
## package 'secr'
## LLsurface.secr.R
## evaluate and plot log likelihood surface for two named beta parameters
## last changed 2010 03 06
############################################################################################

LLsurface.secr <- function (object, betapar = c('g0', 'sigma'), xval = NULL, yval = NULL,
    centre = NULL, realscale = TRUE, plot = TRUE, plotfitted = TRUE, ncores = 1, ...) {

    if (inherits(object, 'list')) {
        temp <- list()
        nsecr <- length(object)
        for (i in 1:nsecr) {
            temp[[i]] <- LLsurface.secr (object[[i]], betapar = betapar, xval=xval, yval=yval,
                centre = centre, realscale = realscale, plot = plot, plotfitted = plotfitted, ...)
        }
        invisible(temp)
    }
    else {
        if (is.null(centre))
             centre <- t(coef(object))['beta',]  ## retains names
        else {
            if (is.null(names(centre)))
                names(centre) <- object$betanames
            else
                if (any(names(centre) != object$betanames))
                    stop ("names of 'centre' do not match 'object$betanames'")
        }
        betaindices <- match(betapar, names(centre))
        if ((length(betapar) != 2) | (any(is.na(betaindices))))
            stop ("requires two named beta parameters")
        if (realscale & any(is.na(match(betapar, names(object$link)))))
            stop ("link function not found - see Notes in help")
        linkx <- ifelse(realscale, object$link[[betapar[1]]], 'identity')
        linky <- ifelse(realscale, object$link[[betapar[2]]], 'identity')
        if (is.null(xval)) {
            betax0 <- centre[betaindices[1]]
            realx0 <- untransform(betax0, linkx)
            xval <- transform (seq(0.8,1.2,0.04) * realx0, linkx)
            xval <- sort(xval)
        }
        else if (realscale) xval <- transform(xval, linkx)
        if (is.null(yval)) {
            betay0 <- centre[betaindices[2]]
            realy0 <- untransform(betay0, linky)
            yval <- transform (seq(0.8,1.2,0.04) * realy0, linky)
            yval <- sort(yval) ## to reverse in case of negative realy0
        }
        else if (realscale) yval <- transform(yval, linky)
        varying <- list(xval,yval)
        names(varying) <- betapar

        grid <- expand.grid(c(as.list(centre[-betaindices]), varying))
        ## restore original order
        grid <- grid[, object$betanames]
        ## drop unnecessary options
        details <- replace(object$details, 'hessian', FALSE)
        details$trace <- FALSE
        details$LLonly <- TRUE

        LL <- function (start) {
            suppressWarnings(
                secr.fit(capthist = object$capthist, model = object$model,
                mask = object$mask, CL = object$CL, detectfn =
                object$detectfn, start = start, link = object$link, fixed
                = object$fixed, timecov = object$timecov, sessioncov =
                object$sessioncov, groups = object$groups, dframe =
                object$dframe, details = details, method =
                object$fit$method, verify = FALSE, ncores = 1)
                )
        }

        cat ('Evaluating log likelihood across grid of', nrow(grid), 'points...\n')
        flush.console()

        if (ncores > 1) {
            clust <- makeCluster(ncores, methods = FALSE, useXDR = .Platform$endian=='big')
            clusterEvalQ(clust, requireNamespace('secr'))   ## faster with this
            temp <- parRapply (clust, grid, LL)
            stopCluster(clust)
        }
        else {
            temp <- apply (grid, 1, LL)
        }

        temp <- matrix(temp, nrow=length(xval))
        if (realscale) {
            xval <- round(untransform(xval, linkx),4)
            yval <- round(untransform(yval, linky),4)
            centre[betapar[1]] <- untransform(centre[betapar[1]], linkx)
            centre[betapar[2]] <- untransform(centre[betapar[2]], linky)
        }

        dimnames(temp) <- list(xval, yval)
        if (plot) {
            contour(x=xval, y=yval, z=temp, xlab=betapar[1], ylab=betapar[2], ...)
            if (plotfitted) {
                points(centre[betapar[1]], centre[betapar[2]], pch = 3)
            }
        }
        invisible(temp)
    }
}

## data(secrdemo)
## LLsurface.secr(secrdemo.0)
