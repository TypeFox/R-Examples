#############################################################################
## package 'secr'
## predict.secr.R
#############################################################################

## 2014-08-19 moved from methods.R
## 2015-09-30 fixpmix now done in secr.lpredictor

############################################################################################
predict.secr <- function (object, newdata = NULL, type = c("response", "link"), se.fit = TRUE,
                          alpha = 0.05, savenew = FALSE, ...) {

    if (is.null(object$fit)) {
        warning ("empty (NULL) object")
        return(NULL)
    }
    type <- match.arg(type)

    if (is.null(newdata)) newdata <- secr.make.newdata (object)

    ## unmashing 2012-07-24
    unmash <- object$details$unmash
    if (object$CL | is.null(unmash)) unmash <- FALSE

    parindices <- object$parindx
    models <- object$model
    smoothsetup <- object$smoothsetup
    if (is.null(smoothsetup)) {
        smoothsetup <- vector(length(models), mode = 'list')
        names(smoothsetup) <- names(models)
    }

    ## drop unused columns 2012-10-24
    vars <- unlist(lapply(models, all.vars))
    ## cover case that h2,h3 not in model 2013-10-28
    mixvar <- switch(object$details$nmix, character(0),'h2','h3')
    usedvars <- c('session', 'g', vars, mixvar)
    newdata <- newdata[,names(newdata) %in% usedvars, drop = FALSE]

    if (object$detectfn %in% c(12,13)) {
        ## experimental parameters not fitted
        ## construct dummies
        parindices$muN <- max(unlist(parindices)) + 1
        parindices$sdN <- max(unlist(parindices)) + 1
        models$muN <- ~1
        models$sdN <- ~1
        object$link$muN <- 'identity'
        object$link$sdN <- 'identity'
    }

    ## allow for fixed beta parameters 2009 10 19, 2014-03-18
    beta <- complete.beta(object)
    beta.vcv <- complete.beta.vcv(object)

    getfield <- function (x) {
        if ((x == 'D') & userD(object)) {
            ## user-supplied density function
            ## return only intercept
            lpred <- matrix(ncol = 2, nrow = nrow(newdata),
               dimnames=list(NULL,c('estimate','se')))
            D0 <- parindices[[x]][1]
            lpred[,1] <- beta[D0]
            lpred[,2] <- beta.vcv[D0,D0]^0.5
            return(lpred)
        }
        else {

            ## smoothsetup argument must be specified if model
            ## includes smooth terms and newdata differs from data
            ## (dframe) used to fit model

            secr.lpredictor (formula = models[[x]], newdata = newdata,
                indx = parindices[[x]], beta = beta, field = x,
                beta.vcv = beta.vcv, smoothsetup = smoothsetup[[x]])
        }
    }
    predict <- sapply (object$realnames, getfield, simplify = FALSE)
  
    z <- abs(qnorm(1-alpha/2))   ## beware confusion with hazard z!
    if (se.fit)  out <- list(nrow(newdata))
    else {
        out <- newdata
        ## add columns for real parameter estimates
        for (varname in object$realnames)
            out[,varname] <- rep(NA,nrow(out))
    }
    
    for (new in 1:nrow(newdata)) {
        lpred  <- sapply (predict, function(x) x[new,'estimate'])
        if (type == "response")
            Xlpred <- Xuntransform(lpred, object$link, object$realnames)

        if (ms(object)) {
            if (!('session' %in% names(newdata))) {
                n.mash <- NULL
            }
            else {
                sess <- newdata[new, 'session']
                n.mash <- attr (object$capthist[[sess]], 'n.mash')
            }
            n.clust <- length(n.mash)
            if (new==1)
                oldnclust <- n.clust
            else if (n.clust != oldnclust)
                warning ("number of mashed clusters varies between sessions")
        }
        else {
            n.mash <- attr (object$capthist, 'n.mash')
            n.clust <- length(n.mash)
        }
        if (unmash)
            n.clust <- 1

        if (se.fit) {
            selpred <- sapply (predict,function(x) x[new,'se'])
            if (type == "response") {
                temp <- data.frame (
                    row.names = object$realnames,
                    link = unlist(object$link[object$realnames]),
                    estimate = Xlpred,
                    SE.estimate = se.Xuntransform (lpred, selpred, object$link, object$realnames),
                    lcl = Xuntransform(lpred-z*selpred, object$link, object$realnames),
                    ucl = Xuntransform(lpred+z*selpred, object$link, object$realnames)
                )
                ## truncate density at zero; adjust for mash(); adjust for telemetry
                if ('D' %in% row.names(temp)) {
                    temp['D', -1][temp['D',-1]<0] <- 0
                    if (!is.null(n.mash)) {
                        temp['D', -1] <- temp['D', -1] / n.clust
                    }
                }
            }
            else {
                temp <- data.frame (
                    row.names = object$realnames,
                    link = unlist(object$link[object$realnames]),
                    estimate = lpred,
                    SE.estimate = selpred,
                    lcl = lpred-z*selpred,
                    ucl = lpred+z*selpred
                )
            }

            # drop non-estimated rows
            if (ms(object)) {
                if ('session' %in% names(newdata))
                    sess <- newdata[new, 'session']
                else
                    sess <- 1
                det <- detector(traps(object$capthist)[[sess]])
            }
            else
                det <- detector(traps(object$capthist))

            ## purge irrelevant real parameters
            if (det == 'telemetry') {
                    rnum <- match(c('D','g0'), row.names(temp))
                    rnum <- rnum[!is.na(rnum)]
                    if (length(rnum)>0)
                        temp <- temp[-rnum,]
            }

            if (nrow(newdata)==1) out <- temp
            else {
                out[[new]] <- temp
                names(out)[new] <- paste (
                        paste(names(newdata),'=', unlist(lapply(newdata[new,],as.character)),
                        sep=' ',collapse=', '),
                    sep=',')
            }
        }
        else { # no SE; terse format
            if (type == "link") {
                out[new, (ncol(newdata)+1) : ncol(out)] <- lpred
            }
            else {
                if ('D' %in% names(Xlpred)) {
                    Xlpred['D'] <- ifelse (Xlpred['D']<0, 0, Xlpred['D'])
                    if (!is.null(n.mash)) {
                        Xlpred['D'] <- Xlpred['D'] / n.clust
                    }
                }
                out[new, (ncol(newdata)+1) : ncol(out)] <- Xlpred
            }
        }
    }
    if (savenew) attr(out, 'newdata') <- newdata
    out
}
############################################################################################

predict.secrlist <- function (object, newdata = NULL, type = c("response","link"), se.fit = TRUE,
                              alpha = 0.05, savenew = FALSE, ...) {
    lapply(object, predict, newdata, type, se.fit, alpha, savenew, ...)
}
############################################################################################
