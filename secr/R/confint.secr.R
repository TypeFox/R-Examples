###############################################################################
## package 'secr'
## confint.secr.R
## last changed 2009 06 11, 2009 07 16 2009 10 20
## 2011 01 31
## 2011 09 28 field argument added to call to secr.lpredictor
## 2011 09 28 optimizer changed to optim from nlm
## 2011 10 21 not for user-defined density
## 2014-08-22 modified to include smoothsetup argument for secr.lpredictor()
## 2014-10-25 updated for NE
## could be sped up by adopting Venzon & Moolgavkar algorithm
## e.g. in Bhat package
###############################################################################

confint.secr <- function (object, parm, level = 0.95, newdata = NULL,
    tracelevel = 1, tol = 0.0001, bounds = NULL, ...) {

    ## profile likelihood interval for estimated secr parameters
    ## cf confint.glm in MASS

    ## interpret character 'parm' as real parameter
    ## interpret numeric 'parm' as beta parameter

    ## case 1 - real parameter not 1:1 beta so require lagrange
    ## case 2 - real parameter but model ~1 so 1:1 beta
    ## case 3 - beta parameter

    #---------------------------------------------------------------------------------------

    profileInterval <- function (parm, ...) {

        predicted <- function (beta) {
            temp <- secr.lpredictor (formula = object$model[[parm]], newdata = newdata,
                indx = object$parindx[[parm]], beta = beta, field = parm,
                smoothsetup = object$smoothsetup[[parm]])[1,'estimate']
            untransform(temp, object$link[[parm]])
        }
        #######################
        ## case 1 - Lagrange

        profileLL.lagrange <- function (gamma, parm) {
            lagrange <- function (beta2, gamma) {
                ## return quantity to be maximized
                templl <- secr.loglikfn (
                    beta       = beta2,
                    parindx    = object$parindx,
                    link       = object$link,
                    fixedpar   = object$fixed,
                    designD    = D.designmatrix,
                    designNE   = NE.designmatrix,
                    design     = object$design,
                    design0    = object$design0,
                    capthist   = object$capthist,
                    mask       = object$mask,
                    detectfn   = object$detectfn,
                    CL         = object$CL,
                    hcov       = object$hcov,
                    groups     = object$groups,
                    details    = details,
                    logmult    = logmult,
                    ncores     = 1,
                    betaw      = max(max(nchar(object$betanames)),8))
                templl - gamma * predicted (beta2)
            }

            ## maximize for fixed gamma (equivalent to fixed 'parm')
            lagrange.fit <- optim (par = object$fit$par, fn = lagrange, gamma = gamma,
                                   hessian = FALSE)
            .localstuff$beta <- lagrange.fit$par
            lp <- - secr.loglikfn (
                beta       = .localstuff$beta,
                parindx    = object$parindx,
                link       = object$link,
                fixedpar   = object$fixed,
                designD    = D.designmatrix,
                designNE   = NE.designmatrix,
                design     = object$design,
                design0    = object$design0,
                capthist   = object$capthist,
                mask       = object$mask,
                detectfn   = object$detectfn,
                CL         = object$CL,
                hcov       = object$hcov,
                groups     = object$groups,
                details    = details,
                logmult    = logmult,
                ncores     = 1,
                betaw      = max(max(nchar(object$betanames)),8))
            ## cat ('gamma ', gamma, '  coef ', lagrange.fit$estimate, '  lp ', lp,
            ##     '  lp - targetLL ', lp-targetLL, '\n')
            lp - targetLL
        }

        #######################
        ## cases 2,3 - fix beta
        
        profileLL3 <- function (x, parm) {
            ## fix required beta parameter to x and evaluate discrepancy
            fb <- object$details$fixedbeta
            if (is.null(fb)) fb <- rep(NA, np)
            fb[parm] <- x
            details$fixedbeta <- fb
            fit <- secr.fit (capthist = object$capthist, mask = object$mask,
                buffer = object$buffer, CL = object$CL, detectfn = object$detectfn,
                start = object$fit$par, binomN= details$binomN,
## link = object$link, fixedpar = object$fixed,  model = object$model,
## timecov = object$timecov, sessioncov = object$sessioncov,
## 2015-12-05
                link = object$link, fixed = object$fixed, model = object$model,
                timecov = object$timecov, sessioncov = object$sessioncov, hcov = object$hcov,
                groups = object$groups, dframe = object$dframe,
                details = details, verify = FALSE, ...)$fit
            - fit$value - targetLL
        }

        #######################
        getlimit <- function (start, step, trialvals) {
            OK <- FALSE
            bound.0 <- start
            f.0 <- profileLL (start, parm)   ## could use fitted LL 31/1/2011

            for (i in trialvals) {
               bound.1 <- start + step * i
               f.1 <- profileLL (bound.1, parm)
               if (prod(sign(c(f.0, f.1)))<0) {  ## different
                   OK <- TRUE
                   break
               }
               else {
                   bound.0 <- bound.1
                   f.0 <- f.1
               }
            }
            if (!OK) {
                if (is.character(parm))
                    stop ("did not find root within ", max(trialvals),
                         " units of MLE")
                else
                    stop ("did not find root within ", max(trialvals),
                         " SE of MLE")
            }
            if (step < 0)
                c(lower = bound.1, upper = bound.0, f.lower = f.1,
                  f.upper = f.0)
            else
                c(lower = bound.0, upper = bound.1, f.lower = f.0,
                  f.upper = f.1)
        }

        #######################

        if (is.numeric(parm)) {
            profileLL <- profileLL3
        }
        else {
            profileLL <- profileLL.lagrange
        }

            estimate <- coef(object, alpha=1-level)[parm,]
            start <- estimate$beta
            se.start <- estimate$SE.beta
            pred <- predict(object)
        if (is.null(bounds)) {
            if (is.numeric(parm))  {
                startlow <- getlimit(estimate$beta, -2 * estimate$SE.beta, c(1,2,4,8))
                startupp <- getlimit(estimate$beta, +2 * estimate$SE.beta, c(1,2,4,8))
            }
            else {
                startlow <- getlimit (0, -1, c(2,5,40,200,1000))
                startupp <- getlimit (0, +1, c(2,5,40,200,1000))
            }
        }
        else {
            startlow <- c(lower = bounds[1],
                          upper = estimate$beta,
                          f.lower = profileLL(bounds[1], parm),
                          f.upper = qchisq(level,1)/2)   ##  1.920729 for level=0.95
            startupp <- c(lower = estimate$beta,
                          upper = bounds[2],
                          f.lower =  qchisq(level,1)/2,  ##  1.920729 for level=0.95
                          f.upper = profileLL(bounds[2], parm))
            validbound <- function(x) prod(sign(x[c('f.lower','f.upper')]))<0
            if (!validbound(startlow))
                stop ("'bounds' does not include lower limit")
            if (!validbound(startupp))
                stop ("'bounds' does not include upper limit")
        }

        temproot <- uniroot (profileLL,
            parm    = parm,
            lower   = startlow['lower'],
            upper   = startlow['upper'],
            f.lower = startlow['f.lower'],
            f.upper = startlow['f.upper'],
            tol     = tol)

        if (is.numeric(parm))
            LCL <- temproot$root
        else
            LCL <- predicted (.localstuff$beta)

        temproot <- uniroot (profileLL,
            parm    = parm,
            lower   = startupp['lower'],
            upper   = startupp['upper'],
            f.lower = startupp['f.lower'],
            f.upper = startupp['f.upper'],
            tol     = tol)

        if (is.numeric(parm))
            UCL <- temproot$root
        else
            UCL <- predicted (.localstuff$beta)

        c(LCL, UCL)
    }
    # end of profileInterval
    #---------------------------------------------------------------------------------------


    memo ('Profile likelihood interval(s)...', tracelevel > 0)

    if (!inherits(object, 'secr'))
        stop ("requires 'secr' object")
    if (userD(object))
        stop ("not implemented for user-defined density function")

    np <- length(object$betanames)  ## number of beta parameters

    ## case 1 - real parameter not 1:1 beta so require lagrange
    ## case 2 - real parameter but model ~1 so 1:1 beta
    ## case 3 - beta parameter

    case <- rep(3, length(parm))

    if (is.character(parm)) {
        OK <- (parm %in% object$realnames)
        if (any(!OK))
            stop ("requested parameter(s) not in model")
        for (i in 1:length(parm))
            if (object$model[[parm[i]]] != ~1)
                case[i] <- 1
            else
                case[i] <- 2
        if (any(case==1)) {
            if (is.null(newdata))
                newdata <- secr.make.newdata (object)[1,, drop = FALSE]  ## default base levels
            if (detector(traps(object$capthist)) %in% c('polygon','polygonX',
                   'transect','transectX','signal','unmarked','presence'))
                logmult <- 0
            else
                logmult <- logmultinom(object$capthist, group.factor(object$capthist,
                    object$groups))

            ## reconstruct density design matrix
            D.modelled <- !object$CL & is.null(object$fixed$D)
            NE.modelled <- is.function(object$details$userdist) & is.null(object$fixed$noneuc)           
            sessionlevels <- session(object$capthist)
            grouplevels <- group.levels(object$capthist, object$groups)
            smoothsetup <- object$smoothsetup
            D.designmatrix <- designmatrix (D.modelled, object$mask, object$model$D,
                                            grouplevels, sessionlevels, object$sessioncov,
                                            smoothsetup$D)
            NE.designmatrix <- designmatrix (NE.modelled, object$mask, object$model$noneuc,
                                            grouplevels, sessionlevels, object$sessioncov,
                                            smoothsetup$noneuc)        
        }
    }
    else {
        if (any ((parm<1) | (parm>np)))
            stop ("invalid beta parameter number")
    }
    #---------------------------------------------------------------------------------------
    targetLL <- - object$fit$value -  qchisq(level,1)/2   # -1.92 for 95% interval
    details <- replace (object$details, 'hessian', FALSE)  ## no need for vcov matrix
    details$trace <- tracelevel > 1
    out <- matrix(nrow = length(parm), ncol = 2)
    if (is.character(parm) & !is.null(bounds)) {
        bounds <- matrix(bounds, ncol = 2)
        for (i in 1:nrow(bounds))
             bounds[i,] <- transform(bounds,object$link[[parm[i]]])
    }

    for (i in 1:length(parm)) {
        parmn <- ifelse (case[i] == 2, match(parm[i], object$betanames), parm[i])
        out[i,] <- profileInterval(parmn)   ## character value if 'real'
        if (case[i] == 2) out[i,] <- untransform(out[i,], object$link[[parm[i]]])
    }
    if (all(case == 3))
        dimnames(out) <- list(object$betanames[parm], c('beta.lcl', 'beta.ucl'))
    else
        dimnames(out) <- list(parm, c('lcl', 'ucl'))
    out
}
############################################################################################
