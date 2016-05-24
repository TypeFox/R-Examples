############################################################################################
## package 'secr'
## esa.R
## last changed 2009 06 29 2009 10 24
## 2010 02 25 dropped details$spherical
## 2010 02 26 fixed for nmix>1
## 2010 03 04 use scaled.detection from functions.R
## 2010 03 09 fixed bug : need to call scaled.detection when !is.null(real)
## 2011-04-04 added noccasions; debugged 2011-04-07
## 2012-10-21 added check to return NA with dettype==13
## 2012-11-13 updated for groups
## 2013-06-06 updated for fixed beta
## 2013-07-20 reparameterize esa,a0
## 2014-08-27 dist2 optional input to integralprwi set to -1
## 2014-09-08 esa adjusted for linearmask cell size
## 2015-10-04 markocc argument for integralprw1
## 2015-11-19 dropped param
############################################################################################

esa <- function (object, sessnum = 1, beta = NULL, real = NULL, noccasions = NULL)

# Return vector of 'a' for given g0, sigma, [z (if hazard fn) ] and session
# detectfn is integer code for detection function 0 = halfnormal, 1 = hazard, 2 = exponential
# 'sessnum' is integer index of session (factor level of the 'session' attribute in capthist)
# object must have at least capthist, mask, detectfn

## strictly doesn't need data, so better to avoid when object not available...
{
    if (ms(object))
        capthists <- object$capthist[[sessnum]]
    else
        capthists <- object$capthist

    if (ms(object$mask))
        mask <- object$mask[[sessnum]]
    else
        mask <- object$mask

    if (is.null(beta) & is.null(real))
        beta <- object$fit$par

    beta <- fullbeta(beta, object$details$fixedbeta)

    trps   <- traps(capthists)  ## need session-specific traps
    if (!(detector(trps) %in% .localstuff$individualdetectors))
        stop ("require individual detector type for esa")
    dettype <- detectorcode(trps)
    n       <- max(nrow(capthists), 1)
    s       <- ncol(capthists)
    constant <- !is.null(noccasions)    ## fix 2011-04-07
    if (is.null(noccasions)) {
        noccasions <- s
    }
    markocc <- markocc(traps(capthists))
    if (is.null(markocc))
        markocc <- rep(1,s)  ## simple marking
    allsighting <- !any(markocc>0)
    if (allsighting) {
        ## drop all zero histories, consider sighting as if marking
        # capthists <- subset(capthists, apply(capthists!=0,1,sum)>0)
        warning("use of allsighting here untested")
        ## markocc[] <- 1
    }

    nmix    <- getnmix (object$details)
    knownclass <- getknownclass(capthists, nmix, object$hcov)

    if (dettype %in% c(3,6)) {
        k <- c(table(polyID(trps)),0)
        K <- length(k)-1
    }
    else if (dettype %in% c(4,7)) {
        k <- c(table(transectID(trps)),0)
        K <- length(k)-1
    }
    else {
        k <- nrow(trps)
        K <- k
    }

    binomN <- object$details$binomN
    m      <- length(mask$x)            ## need session-specific mask...
    cell   <- getcellsize(mask)           ## length or area

    if (constant) {
        ## assume constant
        if (is.null(beta))
            realparval <- detectpar(object)
        else {
            realparval <- makerealparameters (object$design0, beta,
                object$parindx, object$link, object$fixed)  # naive
            realparval <- as.list(realparval)
            names(realparval) <- parnames(object$detectfn)
        }
        a <- cell * sum(pdot(X = mask, traps = trps, detectfn = object$detectfn,
                             detectpar = realparval, noccasions = noccasions))
        return(rep(a,n))
    }
    else {
        if (is.null(beta)) {
            if (is.null(real))
                stop ("requires real parameter values")
            PIA <- rep(1, n * s * K * nmix)    ## nmix added 2010 02 25
            realparval0 <- matrix(rep(real, rep(n,length(real))), nrow = n)   ## UNTRANSFORMED
        }
        else {
            ## allow for old design object
            if (length(dim(object$design0$PIA))==4)
                dim(object$design0$PIA) <- c(dim(object$design0$PIA),1)
            # replaced 2012-09-21
            # PIA <- object$design0$PIA[sessnum,,1:s,,,drop=F]
            PIA <- object$design0$PIA[sessnum,,1:s,1:K,,drop=F]
            ncolPIA <- dim(object$design0$PIA)[2]
            #############################################
            ## trick to allow for changed data 2009 11 20
            ## nmix>1 needs further testing 2010 02 26
            ## NOTE 2010-11-26 THIS LOOKS WEAK
            if (dim(PIA)[2] != n) {
                # 2012-09-21
                PIA <- array(rep(PIA[1,1,,,],n), dim=c(s,K,nmix,n))
                PIA <- aperm(PIA, c(4,1,2,3))   ## n,s,K,nmix
                ncolPIA <- n     ## 2010 02 26
            }
            #############################################

            realparval0 <- makerealparameters (object$design0, beta,
                object$parindx, object$link, object$fixed)  # naive

        }

        ## not compatible with sigmak parameterizations
        Dtemp <- NA
        Xrealparval0 <- reparameterize (realparval0, object$detectfn, object$details,
                                        mask, trps, Dtemp, s)

        ## force to binary 2012-12-17
        ## used <- usage(trps)
        usge <- usage(trps)
        if (is.null(usge)) {
            usge <- matrix(1, nrow = K, ncol = s)
            used <- 1
        }
        else {
            used <- (usge > 1e-10) * 1
        }
        if (any(used == 0))
            PIA <- PIA * rep(rep(t(used),rep(n,s*K)),nmix)
        ncolPIA <- n

        normalize <- object$details$normalize
        if (is.null(normalize))
            normalize <- FALSE
        miscparm <- numeric(4)
        if ((object$detectfn %in% 14:18) & normalize) {
            if (!is.null(object$details$userdist))
                stop("normalization incompatible with userdist")
            miscparm <- c(1,0,1,0)
            if (!is.null(object$details$usecov)) {
                miscparm[2] <- 1
                ## following may fail if some parameters fixed
                alpha2 <- beta[object$parindx[['lambda0']][2]]
                alpha2 <- ifelse (is.na(alpha2),0, alpha2)
                z <- covariates(mask)[,object$details$usecov]
                if (is.null(z))
                    stop("no valid mask covariate")
                mask <- cbind(mask, z * alpha2)
                mask <- as.matrix(mask[,c(1:3,3)]) ## double last col
            }
            else {
                mask <- as.matrix(mask[,c(1:2,2)])
            }
        }
        else
            miscparm[1] <- object$details$cutval
        useD <- FALSE

        ##------------------------------------------
        ## 2014-09-08, 2014-10-17
        if (is.null(object$details$userdist))
            distmat <- -1
        else {

            userdistnames <- getuserdistnames(object$details$userdist)
            if (is.null(covariates(mask)))
                covariates(mask) <- data.frame(row.names = 1:nrow(mask))
            if ('noneuc' %in% userdistnames) {
                covariates(mask)$noneuc <- predictD (object, mask,
                                1, sessnum, parameter = 'noneuc')
            }
            if ('D' %in% userdistnames)
                ## covariates(mask)$D <- D
                stop ("userdist function requiring D not implemented for esa")

            ## pass miscellaneous unmodelled parameter(s) 2015-02-21
            nmiscparm <- length(object$details$miscparm)
            if (nmiscparm > 0) {
                miscindx <- max(unlist(object$parindx)) + (1:nmiscparm)
                attr(mask, 'miscparm') <- coef(object)[miscindx, 1]
            }

            distmat <- valid.userdist (object$details$userdist,
                                       detector(trps),
                                       xy1 = trps,
                                       xy2 = mask,
                                       mask = mask)
        }
        ##------------------------------------------

        ## protection added 2012-10-21
        if (dettype %in% c(13)) {
            return(NA)
        }
        else {
            temp <- .C("integralprw1", PACKAGE = 'secr',
                       as.integer(dettype),
                       as.double(Xrealparval0),
                       as.integer(rep(1,n)),           ## dummy groups 2012-11-13..2013-06-24
                       as.integer(n),
                       as.integer(s),
                       as.integer(k),
                       as.integer(m),
                       as.integer(1),                  ## dummy ngroups 2012-11-13
                       as.integer(nmix),
                       as.integer(knownclass),         ## 2013-04-12
                       as.double(unlist(trps)),
                       as.double(distmat),             ## optional dist2 2014-09-08

                       as.double(usge),
                       as.integer(markocc),
                       as.double(unlist(mask)),
                       as.integer(nrow(Xrealparval0)), ## rows in lookup
                       as.integer(PIA),                ## index of nc*,S,K to rows in realparval0
                       as.integer(ncolPIA),            ## ncol - if CL, ncolPIA = n,
                                                       ## else ncolPIA = 1 or ngrp
                       as.double(cell),                ## cell area (ha) or length (km)
                       as.double(miscparm),
                       as.integer(object$detectfn),
                       as.integer(binomN),             ## 2012-12-18
                       as.integer(useD),
                       a=double(n),
                       resultcode=integer(1)
                       )
            if (temp$resultcode != 0)
                stop ("error in external function 'integralprw1'")
            return(temp$a)
        }
    }
}
############################################################################################
