###############################################################################
## package 'secr'
## mask.check.R
## Evaluate alternative mask buffer width and spacing
## 2010 10 15, 2010-10-17, 2010-10-18, 2010-10-20, 2010-10-31, 2011-06-17
## 2012-11-02 ncores
## 2013-05-03 extended default link for lambda0
## 2015-12-08 uses coefficients of fitted model
###############################################################################

bisect <- function(f, a, b, tol = 1e-6) {
    x1 <- a
    x2 <- b
    f1 <- f(x1)
    f2 <- f(x2)
    while(abs(b - a) > tol) {
       if(f1==f2) {
        b <- x2
        x2 <- x1
        f2 <- f1
        x1 <- b - 0.5 * (b - a)
        f1 <- f(x1)
      } else {
        a <- x1
        x1 <- x2
        f1 <- f2
        x2 <- a + 0.5 * (b - a)
        f2 <- f(x2)
      }
    }
    (a + b) / 2
}

mindistfromedge <- function (mask, traps, maxr=1e5, ntheta = 60, tol=0.0001) {
    ## compare convex hulls
    ## rewritten 2014-08-28; uses inflatechull() in utility.R
    hull <- chull(traps)
    trapshull <- traps[c(hull,hull[1]),]    ## ensure closed
    hull <- chull(mask)
    maskhull <- mask[c(hull,hull[1]),]
    fn <- function (x) {
        all(pointsInPolygon (inflatechull(trapshull, x-tol, ntheta), maskhull))
    }
    maskhull <- inflatechull(maskhull, attr(mask,'spacing')/2, ntheta)
    out <- bisect (fn,0,maxr,tol)
    round(out, round(abs(log10(tol))))
}

mask.check <- function (object, buffers = NULL, spacings = NULL, poly = NULL,
    LLonly = TRUE, realpar = NULL, session = 1, file = NULL, drop = '',
    tracelevel = 0, ncores = 1, ...) {

    dots <- list(...)
    if (!is.null(file)) {
        if (!grepl('.RData',file))
            file <- paste(file, '.RData', sep='')
    }
    refit <- inherits (object,'secr')
    if (!refit & !inherits(object,'capthist'))
        stop ("requires single-session 'capthist' or 'secr' object")
    MS <- ms(object)
    if (MS) {
        if (refit)
            sessnames <-  names(object$capthist)
        else
            sessnames <- names(object)

        ## use alphanumeric session ID
        if (is.numeric(session))
            session <- sessnames[session]

        if (length(session)==1) {
            cat('Session',session, '\n')
        }
    }

    if (MS & (length(session) > 1)) {
        mask.check.outputs <- vector(mode='list')
        if (!LLonly & !is.null(file))
            mask.check.fits <- vector(mode='list')
        if (!is.null(file))
            temp <- paste(tempfile(),'.RData',sep='')
        else
            temp <- NULL
        for (i in session) {
            mask.check.outputs[[i]] <- mask.check (object, buffers, spacings,
                poly, LLonly, realpar, session = i, file = temp, drop,
                tracelevel, ...)
            if (!LLonly & !is.null(file)) {
                load(temp)
                mask.check.fits[[i]] <- mask.check.fit
            }
            if (!is.null(file)) {
                if (LLonly)
                    save (mask.check.outputs, file = file)
                else
                    save (mask.check.outputs, mask.check.fits, file = file)
            }
        }
        mask.check.outputs
    }
    else {                        ## single session
        if (refit) {
            if (!is.call(object$call))
                stop ("requires fitted object from secr 1.5.0 or later with complete call")
            trps <- traps(object$capthist)
            msk <- object$mask
            cpts <- object$capthist
            if (MS) {
                trps <- trps[[session]]
                msk <- msk[[session]]
                cpts <- cpts[[session]]
            }
            newcall <- as.list(object$call)[-1]
            newcall <- replace(newcall, names(dots), dots)
            newcall <- replace(newcall, 'capthist', list(quote(cpts)))
            if (is.null(realpar)) {
                predicted <- predict(object)

                if (MS) {
## 2011-06-28
##                    sessnum <- (1:length(predicted))[grepl(session,
##                        names(predicted))]

## allow for multiple rows; select first
                    sessnum <- (1:length(predicted))[grepl(paste('session =',
                        session), names(predicted))][1]
                    predicted <- predicted[[sessnum]]   ## unreliable?
                }
                realpar <- predicted[,'estimate']
                names(realpar) <- rownames(predicted)
                if (!is.null(newcall$CL))
                    if (newcall$CL == T) {
                        realpar <- realpar[names(realpar) != 'D']
                    }
            }
            one <- is.null(buffers) != is.null(spacings)
            if (is.null(buffers)) {
                distfromedge <- mindistfromedge (msk, trps)
                ## fix 2011-06-17
                if (distfromedge < 1e-10)
                    stop ("specify 'buffers' if detector(s) on edge")
                if (one)
                    buffers <- distfromedge
                else
                    buffers <- c(1, 1.5, 2) * distfromedge
            }
            if (is.null(spacings)) {
                if (one)
                    spacings <- attr(msk, 'spacing')
                else
                    spacings <- c(1, 0.75, 0.5) * attr(msk, 'spacing')
            }

        }
        else {
            if (is.null(buffers))
                stop ("requires 'buffers'")
            if (is.null(spacings))
                stop ("requires 'spacings'")
            cpts <- object
            if (MS) {
                cpts <- cpts[[session]]
            }
            if (tracelevel>0) verify(cpts)
            newcall <- list(capthist = quote(cpts))
            newcall <- replace(newcall, names(dots), dots)
            trps <- traps(cpts)
        }
        newcall <- replace(newcall, 'trace', tracelevel>1)
        newcall <- replace(newcall, 'verify', FALSE)
        nbuffer <- length(buffers)
        nspacing <- length(spacings)

        if (LLonly) {
            defaultlink <- list(D='log', g0='logit', lambda0='log', sigma='log', z='log',
                w='log', pID='logit', beta0='identity', beta1='neglog',
                sdS='log', b0='log', b1='neglog', pmix='logit')
            if (detector(trps) %in% .localstuff$countdetectors) defaultlink$g0 <- 'log'
            newcall$link <- defaultlink
            
            ## use betas from fitted model 2015-12-08
            if (inherits(object, 'secr')) {
                startbeta <- coef(object)$beta
            }
            else {
                if (is.null(realpar))
                    stop ("'LLonly' requires either 'realpar' or fitted model")
                startbeta <- Xtransform(unlist(realpar), defaultlink, names(realpar))
            }
            newcall <- replace(newcall, 'start', list(start=startbeta))
            
            if (is.null(newcall$details)) {
                newcall$details <- list(LLonly = TRUE)
            }
            else
                newcall$details <- replace(newcall$details, 'LLonly', TRUE)
            mask.check.output <- array(dim=c(nbuffer, nspacing))
            dimnames(mask.check.output) <- list(buffer=buffers, spacing=spacings)
            cat('Computing log likelihoods... \n')
        }
        else {
            mask.check.output <- array(dim=c(nbuffer, 4, nspacing))
            dimnames(mask.check.output) <- list(buffer=buffers,
                c('esa','LL','D-hat','se(D-hat)'), spacing=spacings)
            if (!is.null(file)) {
                mask.check.fit <- vector(mode = 'list')
                class(mask.check.fit) <- c('list','secrlist')
            }
            cat('Fitting models... \n')
        }
        if (!is.null(newcall$model))
            if (any(c('session','Session') %in% as.character(newcall$model)) |
                !is.null(newcall$sessioncov)) {
                stop ("mask.check cannot handle models with session structure")
            }

        if (ncores > 1) {
            sfit <- function (ij) {
                spacing <- ij[1]
                buffer <- ij[2]
                print(spacing)
                msk <- make.mask(traps = trps, buffer = buffer, spacing = spacing,
                                 type = 'trapbuffer', poly = poly)
                newcall <- replace(newcall, 'mask', list(quote(msk)))
                if (tracelevel<2) old <- options(warn=-1)
                if (LLonly) {
                    do.call(secr.fit, newcall)
                }
                else {
                    newfit <- do.call(secr.fit, newcall)
                    if (newfit$CL) {
                        tempij <- as.matrix(derived(newfit))
                        c(tempij[1,1], -newfit$fit$minimum, tempij[2,1:2])
                    }
                    else {
                        c(esa(newfit)[1], -newfit$fit$minimum,
                          unlist(predict(newfit)['D', c('estimate','SE.estimate')]))
                    }
                }
            }
            ij.df <- expand.grid (i=spacings, j=buffers)
            clust <- makeCluster(ncores, methods = FALSE, useXDR = .Platform$endian=='big')
            clusterEvalQ(clust, requireNamespace('secr'))
            output <- parRapply(clust, ij.df, sfit)
            stopCluster(clust)

            if (LLonly)
                mask.check.output[] <- t(matrix(output, nrow=nspacing, ncol=nbuffer))
            else
                mask.check.output[] <- aperm(array(output, dim=c(4, nspacing, nbuffer)),c(3,1,2))
        }
        else {
            ## old code for ncores = 1
            flush.console()
            for (j in 1:nspacing) {
                spacing <- spacings[j]
                for (i in 1:nbuffer) {
                    buffer <- buffers[i]
                    msk <- make.mask(trps, buffer, spacing, type = 'trapbuffer',
                                     poly = poly)
                    ## quote passes name... avoids bulky call
                    newcall <- replace(newcall, 'mask', list(quote(msk)))
           
                    if (tracelevel<2) old <- options(warn=-1)
                    if (LLonly) {
                        mask.check.output[i,j] <- do.call(secr.fit, newcall)
                    }
                    else {
                        newfit <- do.call(secr.fit, newcall)
                        if (newfit$CL) {
                            temp <- as.matrix(derived(newfit))
                            mask.check.output[i,1,j] <- temp[1,1]
                            mask.check.output[i,2,j] <- -newfit$fit$minimum
                            mask.check.output[i,3:4,j] <- temp[2,1:2]
                        }
                        else {
                            mask.check.output[i,1,j] <- esa(newfit)[1]
                            mask.check.output[i,2,j] <- -newfit$fit$minimum
                            mask.check.output[i,3:4,j] <- unlist(predict(newfit)
                                                                 ['D', c('estimate','SE.estimate')])
                        }
                    }
                    options(old)
                if (tracelevel>0)
                    cat('Completed buffer', buffers[i], 'spacing',
                        spacings[j], date(), '\n')
                    flush.console()
                    if (!is.null(file)) {
                        if (LLonly) {
                            save (mask.check.output, file = file)
                        }
                        else {
                            modelname <- paste('buffer', buffers[i], 'spacing',
                                               spacings[j], sep = '')
                            newfit <- trim(newfit, drop)
                            mask.check.fit[[modelname]] <- newfit
                            save (mask.check.output, mask.check.fit, file = file)
                        }
                    }
                }
            }
        }  ## end of old code for ncores=1

        if (!LLonly)
        {
            if ((nbuffer == 1) & (nspacing>1)) ## drop first dim
                mask.check.output <- t(mask.check.output[1,,])
            else
                if (nspacing == 1)             ## drop third dim
                    mask.check.output <- mask.check.output[,,1]
        }
        mask.check.output
    }
}
###############################################################################
