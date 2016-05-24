############################################################################################
## package 'secr'
## regionN.R
## population size in an arbitrary region
## 2011-08-18 (fixed expected.n)
## 2011-09-26 fixed multi-session bug in region.N
## 2011-10-19 adjustments for speed, observe se.N
## 2011-10-19 slowness is due to call of integralprw1 in sumDpdot, esp in betaRN
## 2011-10-20 minor editing
## 2011-10-21 predictD moved to Dsurface.R
## 2012-04-18 bug fixed: nlowerbound and RN.method ignored when nsess>1
## 2012-05-13 added explicit 'poisson' option in region.N for computation of RN

## 2012-08-07 potentially extend to groups by looping over sessions _and_ groups
## 2014-03-19 pooled.RN
## 2014-08-27 dist2 optional input to integralprwi set to -1
## 2014-09-08 renamed cellarea to cell, inclusive of cell length for linearmask
## 2014-09-08 renamed regionarea to regionsize
## 2014-09-08 updated with reparameterize
## 2014-09-08 prepared for secrlinear linearmask models, not activated
## 2015-10-04 markocc argument for integralprw1
## 2015-11-19 dropped param argument for integralprw1
############################################################################################

region.N <- function (object, region = NULL, spacing = NULL, session = NULL,
    group = NULL, se.N = TRUE, alpha = 0.05, loginterval = TRUE,
    keep.region = FALSE, nlowerbound = TRUE, RN.method = 'poisson',
    pooled.RN = FALSE) {

    ## Notes
    ## se.N = FALSE returns scalar N

    ###########################################################
    ## for gradient of E.N wrt density betas
    betaEN <- function (betaD, object, regionmask, group, session) {
        ## regionmask is a mask (no need for spacing)
        ## assume single session
        ## indx identifies beta parameters for density D
        object$fit$par[indx] <- betaD
        region.N(object, regionmask, spacing = NULL, session = session,
            group = group, se.N = FALSE, keep.region = FALSE,
            pooled.RN = FALSE)
    }
    ###########################################################
    ## for gradient of R.N wrt all betas
    betaRN <- function (beta, object, regionmask) {
        ## regionmask is a mask (no need for spacing)
        ## assume single session
        ## n, cell, sessnum global
        object$fit$par <- beta
        D <- predictD(object, regionmask, group, session, parameter = 'D')
        noneuc <- predictD(object, regionmask, group, session, parameter = 'noneuc')
        n + sumDpdot (object, sessnum, regionmask, D, noneuc, cell,
                  constant = FALSE, oneminus = TRUE, pooled = pooled.RN)[1]
    }
    ###########################################################

    if (is.null(region)) {
        region <- object$mask
        ## warning ("using entire mask as region")
    }

    if (!all(session %in% session(object$capthist)))
        stop ("session incompatible with object ")

    if (is.null(group))
        group <- 1

    if (is.null(session))
        session <- session(object$capthist)

    if (!ms(object))
        pooled.RN <- FALSE

    if (pooled.RN)
        ## use only model for first session
        session <- session[1]

    ####################################################################
    ## if N requested for multiple sessions,
    ## call region.N recursively for each session
    nsess <- length(session)
    if (nsess > 1) {
        ## predict for each session
        out <- vector('list')
        for (sess in session) {
            if (ms(region))
                tempregion <- region[[sess]]
            else
                tempregion <- region
            out[[sess]] <- region.N (object, region = tempregion,
                spacing = spacing, session = sess, group = group,
                se.N = se.N, alpha = alpha, loginterval = loginterval,
                keep.region = keep.region, nlowerbound = nlowerbound,
                RN.method = RN.method, pooled.RN = FALSE)
        }
        out
    }

    ####################################################################
    ## otherwise, this is a non-recursive call for one session...
    else {
        if (ms(object$mask))
            mask <- object$mask[[session]]
        else
            mask <- object$mask
        masktype <- if (inherits(mask, "linearmask")) "linearmask" else "mask"

        ########################################################
        ## if necessary, convert vector region to raster
        if (inherits(region, 'mask')) {
            ## includes linearmask
            if (ms(region))
                 regionmask <- region[[session]]
            else
                regionmask <- region
        }
        else {
            if (is.null(spacing)) {
                ## use mask spacing by default
                spacing <- spacing(mask)
            }
            if ((inherits(region, 'SpatialPolygonsDataFrame') & (masktype == "mask")) |
                (inherits(region, 'SpatialLinesDataFrame') & (masktype == "linearmask"))
                ) {
                bbox <- bbox(region)
            }
            else {
                bbox <- apply(region, 2, range)
            }
            if (masktype == "mask") {
                regionmask <- make.mask(bbox, poly = region, buffer = 0,
                                        spacing = spacing, type = 'polygon',
                                        check.poly = FALSE)
            }
            else {
                ## dodge circular dependence and anticipate release of secrlinear
                ## 2014-09-09
                stop("construction of linear masks not yet supported in region.N")
                ## if (!requireNamespace("secrlinear", quietly=TRUE))
                ##     stop ("could not load secrlinear")
                ## regionmask <- secrlinear::read.linearmask(data = region,
                ##     spacing = spacing, spacingfactor = attr(mask, "spacingfactor"))
            }
        }

        ## 2015-01-16
        if (is.matrix(object$details$userdist))
            if (ncol(object$details$userdist) != nrow(regionmask))
                stop("userdist matrix incompatible with region mask")

        #######################################################

        ## region now inherits from mask, so has area attribute
        if (inherits(regionmask, "linearmask"))
            cell <- attr(regionmask, 'spacing') / 1000    ## km
        else
            cell <- attr(regionmask, 'area')              ## ha
        regionsize <- nrow(regionmask) * cell

        if (ms(object)) {
            if (pooled.RN)
                n <- sum(sapply(object$capthist, nrow))
            else
                n <- nrow(object$capthist[[session]])
        }
        else
            n <- nrow(object$capthist)
        sessnum <- match (session, session(object$capthist))

        #######################################################
        ## for conditional likelihood fit,
        if (object$CL) {
            temp <- derived(object, se.D = se.N) ## inefficient as repeats for each sess
            if (!is.data.frame(temp))
                temp <- temp[[session]]
            D <- temp['D', 'estimate']
            seD <- temp['D', 'SE.estimate']
            EN <- D * regionsize
            if (!se.N) return (EN)    ## and stop here
            seEN <- seD * regionsize
        }

        #######################################################
        ## for full likelihood fit...
        else {
            if (is.null(object$model$D) | is.null(object$link$D))
                stop ("model or link function not found in object")

            if ((object$model$D == ~1) & !userD(object)) {
                predicted <- predict(object)
                if (!is.data.frame(predicted))
                    predicted <- predicted[[1]]
                D <- predicted['D','estimate']
                seD <- predicted['D', 'SE.estimate']

                EN <- D * regionsize
                if (!se.N) return (EN)    ## and stop here
                seEN <- seD * regionsize
            }
            else {
                D <- predictD (object, regionmask, group, session, parameter = 'D')
                EN <- sum(D) * cell
                if (!se.N) return (EN)    ## and stop here
                indx <- object$parindx$D
## simple gradient failed in test 2011-10-20
#               dENdphi <- gradient (object$fit$par[indx],
#                    betaEN, object = object, region = region, session =
#                    session, group = group)
                dENdphi <- nlme::fdHess (object$fit$par[indx],
                    betaEN, object = object, region = regionmask, group = group,
                    session = session)$gradient
                beta.vcv <- object$beta.vcv[indx,indx]
                seEN <- (dENdphi %*% beta.vcv %*% dENdphi)^0.5
            }
        }

        #######################################################################
        ## realised N
        ## only makes sense for individual detectors (not unmarked or presence)
        ## assume if we have got this far that SE is required
        ## amended 2011-09-26
        if (ms(object))
            det <- detector(traps(object$capthist)[[session]])
        else
            det <- detector(traps(object$capthist))

        if (det %in% .localstuff$individualdetectors) {
            noneuc <- predictD (object, regionmask, group, session, parameter = 'noneuc')
            RN.method <- tolower(RN.method)
            if (RN.method == 'mspe') {
                notdetected <- sumDpdot (object, sessnum, regionmask, D, noneuc,
                    cell, constant = FALSE, oneminus = TRUE, pooled = pooled.RN)[1]
                RN <- n + notdetected
                ## evaluate gradient of RN wrt betas at MLE
                dNdbeta <- nlme::fdHess (object$fit$par, betaRN, object = object,
                    region = regionmask)$gradient
                ## compute variance from gradient & vcv
                pdotvar <- dNdbeta %*% object$beta.vcv %*% dNdbeta
                seRN <- (notdetected + pdotvar)^0.5
            }
            else if (RN.method == 'poisson') {
                notdetected <- sumDpdot (object, sessnum, regionmask, D, noneuc,
                    cell, constant = FALSE, oneminus = TRUE, pooled = pooled.RN)[1]
                RN <- n + notdetected
                seRN <- (seEN^2 - EN)^0.5
            }
            ## RN.method = 'EN'
            else {
                RN <- EN
                seRN <- (seEN^2 - EN)^0.5
            }
        }
        else { RN <- NA; seRN <- NA }

# suppress 2011-11-10
#            ## additional flourish - compute expected n
#            En <- sumDpdot (object, sessnum, regionmask, D, attr(regionmask,'area'),
#                 constant = FALSE, oneminus = FALSE)[1]
#        else { RN <- NA; seRN <- NA; En <- NA }
        #######################################################################

        temp <- data.frame(
            row.names = c('E.N','R.N'),
            estimate = c(EN,RN),
            SE.estimate = c(seEN,seRN))
        ## lower bound added 2011-07-15
        ## vector (0,n) means apply to R.N not E.N
        if (nlowerbound)
            temp <- add.cl (temp, alpha, loginterval, c(0, n))
        else
            temp <- add.cl (temp, alpha, loginterval, c(0, 0))
        temp$n <- rep(n, nrow(temp))
 #       temp$E.n <- rep(round(En,2), nrow(temp))
        ## 2014-11-12
        if (inherits(region, 'linearmask'))
            attr(temp, 'regionsize') <- masklength(region)
        else
            attr(temp, 'regionsize') <- maskarea(region) ## nrow(region) * attr(region, 'area')
        if (keep.region)
            attr(temp, 'region') <- region
        temp
    }
}

############################################################################################
## 2011-05-05
############################################################################################

## modelled on esa.R

sumDpdot <- function (object, sessnum = 1, mask, D, noneuc, cell, constant = TRUE,
                      oneminus = FALSE, pooled = FALSE)

# Return integral for given model and new mask, D
# 'sessnum' is integer index of session (factor level of the 'session' attribute in capthist)
# object must have at least capthist, detectfn
# D should be scalar or vector of length nrow(mask)

# if 'constant' a much simplified calculation is used, assuming
# constant detection and density, and full detector usage

{
    if (ms(object))
        capthists <- object$capthist[[sessnum]]
    else
        capthists <- object$capthist

    if (ms(mask))
        mask <- mask[[sessnum]]

    ## allow for fixed beta parameters 2014-03-18
    beta <- complete.beta(object)

    n       <- max(nrow(capthists), 1)
    s       <- ncol(capthists)
    noccasions <- s

    ## 2014-03-19 pooling option
    if (pooled & ms(object))
        trps <- do.call(rbind, c(traps(object$capthist), list(addusage = TRUE)))
    else
        trps   <- traps(capthists)  ## use session-specific traps
    if (!(detector(trps) %in% .localstuff$individualdetectors))
        stop ("require individual detector type for sumDpdot")

    dettype <- detectorcode(trps)
    nmix    <- getnmix(object$details)
    knownclass <- getknownclass(capthists, nmix, object$hcov)

    ##############################################
    ## marking occasions 2015-10-04
    markocc <- markocc(traps(capthists))
    if (is.null(markocc))
        markocc <- rep(1,s)
    ##############################################

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
    if (constant) {
        if (is.null(beta))
            real <- detectpar(object)
        else {

            real <- makerealparameters (object$design0, beta,
                object$parindx, object$link, object$fixed)  # naive
            real <- as.list(real)
            names(real) <- parnames(object$detectfn)
        }
        a <- cell * sum(pdot(X = mask, traps = trps, detectfn = object$detectfn,
                             detectpar = real, noccasions = noccasions))
        return(a * D)
    }
    else {

        ## allow for old design object
        if (length(dim(object$design0$PIA))==4)
            dim(object$design0$PIA) <- c(dim(object$design0$PIA),1)
        PIA <- object$design0$PIA[sessnum,,1:s,,,drop=F]
        ncolPIA <- dim(object$design0$PIA)[2]
        #############################################
        ## trick to allow for changed data 2009 11 20
        ## nmix>1 needs further testing 2010 02 26
        ## NOTE 2010-11-26 THIS LOOKS WEAK
        if (dim(PIA)[2] != n) {
            PIA <- array(rep(PIA[1,1,,,],n), dim=c(s,K,nmix,n))
            PIA <- aperm(PIA, c(4,1,2,3))   ## n,s,K,nmix
            ncolPIA <- n     ## 2010 02 26
        }
        #############################################

        realparval0 <- makerealparameters (object$design0, beta,
            object$parindx, object$link, object$fixed)  # naive

        ## 2014-09-08
        if (!is.null(object$fixed$D))
            Dtemp <- object$fixed$D
        else if (object$CL)
            Dtemp <- NA
        else
            Dtemp <- D[1]
                
        Xrealparval0 <- reparameterize (realparval0, object$detectfn, object$details,
                                        mask, trps, Dtemp, s)

        usge <- usage(trps)
        if (is.null(usge)) {
            usge <- matrix(1, nrow = K, ncol = s)
            used <- 1
        }
        else {
            used <- (usge > 1e-10) * 1
        }
        if (any(used==0))
        PIA <- PIA * rep(rep(t(used),rep(n,s*K)),nmix)
        ncolPIA <- n

        miscparm <- numeric(4)
        miscparm[1] <- object$details$cutval

        ## add density as third column of mask
        if (!(length(D) %in% c(1,nrow(mask))))
            stop ("D does not match mask in sumDpdot")

        ##------------------------------------------
        ## 2014-09-08
        if (is.null(object$details$userdist))
            distmat <- -1
        else {

            userdistnames <- getuserdistnames(object$details$userdist)
            if (is.null(covariates(mask)))
                covariates(mask) <- data.frame(row.names = 1:nrow(mask))
            if ('noneuc' %in% userdistnames)
                covariates(mask)$noneuc <- noneuc
            if ('D' %in% userdistnames)
                covariates(mask)$D <- D

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

        if (length(D) == 1) {
            useD <- FALSE
        }
        else {
            mask <- cbind (mask, D)
            useD <- TRUE
        }

        temp <- .C("integralprw1", PACKAGE = 'secr',
            as.integer(dettype),
            as.double(Xrealparval0),
            as.integer(rep(1,n)),           ## dummy groups 2012-11-13; 2013-06-24
            as.integer(n),
            as.integer(s),
            as.integer(k),
            as.integer(m),
            as.integer(1),                  ## dummy ngroups 2012-11-13
            as.integer(nmix),
            as.integer(knownclass),
            as.double(unlist(trps)),
            as.double(distmat),             ## optional dist2 2014-09-08
            as.double(usge),
            as.integer(markocc),            ## 2015-10-04
            as.double(as.numeric(unlist(mask))),
            as.integer(nrow(Xrealparval0)), ## rows in lookup
            as.integer(PIA),                ## index of nc*,S,K to rows in realparval0
            as.integer(ncolPIA),            ## ncol - if CL, ncolPIA = n, else ncolPIA = 1 or ngrp
            as.double(cell),
            as.double(miscparm),
            as.integer(object$detectfn),
            as.integer(binomN),
            as.integer(useD),
            a = double(n),
            resultcode = integer(1)
       )
       if (temp$resultcode != 0)
           stop ("error in external function 'integralprw1'")

       ## constant density case, D not passed to integralprw1
       if (length(D) == 1) {
           temp$a <- temp$a * D
       }

       if (oneminus) {
           sumD <- ifelse (length(D) == 1, D * nrow(mask), sum(D))
           return(sumD * cell - temp$a)
       }
       else
           return(temp$a)
    }
}
############################################################################################


expected.n <- function (object, session = NULL, group = NULL, bycluster = FALSE,
                        splitmask = FALSE) {

    ## Note
    ## splitmask toggles between two methods for clustered detectors:
    ## 1. is to integrate over the whole mask, restricting detection to each cluster in turn
    ## 2. is to split mask into Dirichlet subregions by distance to detector centre
    ## and to integrate over all detectors, assuming those far away will never detect from
    ## within a subregion to which they do not belong
    ## Probably, 1. is more robust

    if (!all(session %in% session(object$capthist)))
        stop ("session incompatible with object")

    if (!is.null(group))
        stop ("not yet working for groups")

    if (!all(group %in% interaction (object$groups)))
        stop ("unrecognised groups")

    if (is.null(session))
        session <- session(object$capthist)

    ####################################################################
    ## if En requested for multiple sessions,
    ## call En recursively for each session
    nsess <- length(session)
    if (nsess > 1) {
        ## predict for each session
        out <- vector('list')
        for (sess in session) {
            out[[sess]] <- expected.n (object, sess, group, bycluster, splitmask)
        }
        out
    }

    ####################################################################
    ## otherwise, this is a non-recursive call for one session...
    else {

        if (ms(object$mask))
            mask <- object$mask[[session]]
        else
            mask <- object$mask
        if (inherits(mask, "linearmask"))
            cell <- attr(mask, 'spacing') / 1000    ## km
        else
            cell <- attr(mask, 'area')              ## ha
        if (ms(object)) {
            n <- nrow(object$capthist[[session]])
            trps <- traps(object$capthist[[session]])
        }
        else {
            n <- nrow(object$capthist)
            trps <- traps(object$capthist)
        }
        sessnum <- match (session, session(object$capthist))

        #######################################################
        ## for conditional likelihood fit,
        if (object$CL) {
            temp <- derived(object, se.D = FALSE) ## inefficient as repeats for each sess
            if (!is.data.frame(temp))
                temp <- temp[[session]]
            D <- temp['D', 'estimate']
        }
        #######################################################
        ## for full likelihood fit...
        else {
            if (is.null(object$model$D) | is.null(object$link$D))
                stop ("model or link function not found in object")

            if (object$model$D == ~1) {
                predicted <- predict(object)
                if (!is.data.frame(predicted))
                    predicted <- predicted[[1]]
                D <- rep(predicted['D','estimate'], nrow(mask))
            }
            else {
                D <- predictD (object, mask, group, session, parameter = 'D')
            }
        }
        #######################################################

        if (is.function(object$details$userdist)) {
          if (object$model$noneuc == ~1) {
            predicted <- predict(object)
            if (!is.data.frame(predicted))
              predicted <- predicted[[1]]
            noneuc <- rep(predicted['noneuc','estimate'], nrow(mask))
          }
          else {
            noneuc <- predictD (object, mask, group, session, parameter = 'noneuc')
          }
        }
        else noneuc <- rep(NA, nrow(mask))

        #################################################################
        if (bycluster) {
            centres <- cluster.centres(trps)
            nclust <- nrow(centres)
            out <- numeric (nclust)
            if (is.null(attr(trps, 'cluster'))) {
                clusterID(trps) <- 1:nclust
            }
            if (splitmask) {
                cluster <- nearesttrap (mask, centres)
                mask <- split (mask, cluster)
                D <- split(D, cluster)
                noneuc <- split(noneuc, cluster)
            }
            for (i in 1:nclust) {
                if (splitmask) {
                    out[i] <- sumDpdot(object = object, sessnum = sessnum,
                        mask=mask[[i]], D = D[[i]], noneuc = noneuc[[i]], cell = cell,
                        constant = FALSE, oneminus = FALSE)[1]
                }
                else {
                    temptrap <- subset(trps, subset = as.numeric(clusterID(trps)) == i)
                    if (ms(object))
                        traps(object$capthist[[sessnum]]) <- temptrap
                    else
                        traps(object$capthist) <- temptrap

                    out[i] <- sumDpdot(object = object, sessnum = sessnum,
                        mask=mask, D = D, noneuc = noneuc, cell = cell,
                        constant = FALSE, oneminus = FALSE)[1]
                }
            }
            out
        }
        else {
            sumDpdot (object, sessnum, mask, D, noneuc, cell,
             constant = FALSE, oneminus = FALSE)[1]
        }
        #################################################################
    }
}

