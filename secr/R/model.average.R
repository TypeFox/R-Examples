############################################################################################
## package 'secr'
## model.average.R
## Last changed 2010 02 17 fixed use of newdata (previously ignored)
## 2010 03 09 fixed sort order of collate with sessioncov
## 2011-11-05 covar = TRUE bugfix
## 2013-06-10 criterion
## 2013-06-10 MATA
## 2015-09-30 added annotations
## 2015-09-30 fixpmix now done in secr.lpredictor
############################################################################################

MATA <- function (wt, est, se, alpha) {
    ## Turek and Fletcher 2012 "model-averaged tail area" interval
    fn <- function (theta, lcl = TRUE) {
        sum(wt * pnorm((est-theta) / se, lower.tail = lcl )) - alpha/2
    }
    lcl <- uniroot (fn, interval = c(min(est - 10*se), max(est)), lcl = FALSE)$root
    ucl <- uniroot (fn, interval = c(min(est), upper = max(est + 10*se)), lcl = TRUE)$root
    cbind(lcl, ucl)
}

model.average <- function (..., realnames = NULL, betanames = NULL,
    newdata = NULL, alpha = 0.05, dmax = 10, covar = FALSE, average =
    c('link', 'real'), criterion = c('AICc','AIC'), CImethod =
    c('Wald', 'MATA')) {

    ########
    ## SETUP

    ## match character arguments
    CImethod <- match.arg(CImethod)
    criterion <- match.arg(criterion)
    average <- match.arg(average)

    ## coerce input models to secrlist
    object <- secrlist(...)
    if (!is.null(names(object)))
        modelnames <- names(object)
    else
        modelnames <- as.character(match.call(expand.dots=FALSE)$...)
    
    ## checks
    if ( any (!sapply(object, function (x) inherits(x, 'secr'))) )
        stop ("require fitted 'secr' objects")
    if ( length(object) < 2 )
        warning ("only one model")
    if (!is.list(object) | !inherits(object[[1]], 'secr'))
        stop ("object must be secr or list of secr")

    ## preliminaries
    type <- 'real'                     ## default
    parnames <- object[[1]]$realnames  ## default
    links <- object[[1]]$link
    if (!is.null(realnames))
        parnames <- realnames
    else if (!is.null(betanames)) {
        type <- 'beta'
        average <- 'beta'   ## override
        parnames <- betanames
    }
    np <- length(parnames)
    nsecr <- length(object)
    
    ## rudimentary checks for compatible models
    if (nsecr > 1) {
        objnames <- function(i) switch (type, real=object[[i]]$realnames, beta=object[[i]]$betanames)
        test <- sapply (2:nsecr, function(i) sum(match(parnames, objnames(i), nomatch=0)>0) == np)
        if (!all(test))
            stop ("requested parameters not found in all models, ",
                  "or models incompatible")
    }

    ## NEWDATA construct if not provided
    if (is.null(newdata)) {
        #############################################################
        ## form unified 'newdata' containing all necessary predictors
        tempnewdata <- lapply (object, secr.make.newdata)
        ## extract list of columns from all components of 'newdata'
        ## modified 2010 02 14
        column.list <- list(0)
        for (i in 1:nsecr) column.list[[i]] <- as.list(tempnewdata[[i]])
        column.list <- unlist(column.list, recursive = FALSE)

        column.list <- column.list[unique(names(column.list))]
        column.list <- lapply(column.list, unique)
        common <- names(column.list)[names(column.list) %in% names(newdata)]
        column.list[common] <- newdata[common]   ## user input

        ## modified 2010 02 23 to match levels of session covariates to session
        ## not tested with different session covariates in diff sessions
        sessioncovs <- lapply(object, function(x)
            if(!is.null(x$sessioncov)) data.frame(session=session(x$capthist), x$sessioncov)
            else NULL)
        sessioncovs <- sessioncovs[!sapply(sessioncovs, is.null)]
        scn <- as.vector(sapply(sessioncovs, names))
        scn <- match(unique(scn),scn)
        sessioncovs <- as.data.frame(sessioncovs)[,scn]
        sessioncovs <- as.data.frame(sessioncovs[!sapply(sessioncovs, is.null)])
        sessioncovnames <- unlist(lapply(object, function(x) names(x$sessioncov)))
        sessioncovariate <- names(column.list) %in% sessioncovnames
        newdata <- expand.grid (column.list[!sessioncovariate])
        ## replaced 2010 03 09 to retain sort order
        #   if (nrow(sessioncovs)>0) newdata <- merge(newdata, sessioncovs)
        if (nrow(sessioncovs)>0) {
            for (i in names(sessioncovs)) {
                if (i != 'session') newdata[,i] <- sessioncovs[newdata$session,i]
            }
        }
    }
    
    nr <- nrow(newdata)
    rownames <- apply(newdata, 1, function(x) paste(names(newdata), '=', x, sep='', collapse=','))

    ## VECTOR OF MODEL WEIGHTS

    AICdf <- AIC(object, sort = FALSE, dmax = dmax, criterion = criterion)
    wt <- unlist(AICdf[,paste(criterion, 'wt', sep = '')])

    ## AVERAGE MODELS

    if (type == 'beta') {
        nr   <- 1
        ests <- lapply (object, function(x) coef(x)[parnames, 'beta'])
        ests <- array(unlist(ests), dim=c(nr, np, nsecr))
        estse <- lapply (object, function(x) coef(x)[parnames, 'SE.beta'])   ## for MATA
        estse <- array(unlist(estse), dim=c(nr, np, nsecr))
        ests <- array(unlist(ests), dim=c(nr, np, nsecr))
        estse <- array(unlist(estse), dim=c(nr, np, nsecr))
        wtd <- sweep(ests, MARGIN = 3, STATS = wt, FUN = '*')
        wtd <- apply(wtd, 1:2, sum)
        ## fails if no dimnames on x$beta.vcv:
        vcvs <- lapply (object, function(x) { vcov(x)[parnames, parnames] })
        vcv1 <- array(unlist(vcvs), dim=c(np, np, nsecr))  ## between beta parameters
        vcv <- sweep (vcv1, MARGIN = 3, STATS = wt, FUN = '*')
        vcv <- apply(vcv, 1:2, sum)
        sewtd  <- diag(vcv)^0.5
    }
    else { ## type == 'real'
        getLP <- function (object1) {  ## predicted values of real parameters
            getfield <- function (x) {
                 secr.lpredictor (formula = object1$model[[x]], newdata = newdata,
                    indx = object1$parindx[[x]], beta = object1$fit$par,
                    beta.vcv = object1$beta.vcv, field = x,
                    smoothsetup = object1$smoothsetup[[x]])
            }
            sapply (names(object1$model), getfield, simplify = FALSE)
        }
        predicted <- lapply (object, getLP)

        ests <- lapply (predicted, function(x) lapply(x[parnames], function(y) y[, 'estimate'] ))
        estse <- lapply (predicted, function(x) lapply(x[parnames], function(y) y[, 'se'] ))
        if (average == 'real') {
            ests <- lapply (ests, function (x) Xuntransform(unlist(x),
                varnames=rep(parnames, rep(nr, np)), linkfn=links))
            estse <- mapply (function (x,se.x) se.Xuntransform(unlist(x), unlist(se.x),
                varnames=rep(parnames, rep(nr, np)), linkfn=links), ests, estse)
            vcvs <- lapply (object, vcov, realnames = parnames, newdata=newdata)
        }
        else {
            vcvs <- lapply (predicted, function(x) lapply(x[parnames], function(y) attr(y, 'vcv')))
        }
        ests <- array(unlist(ests), dim=c(nr, np, nsecr))
        estse <- array(unlist(estse), dim=c(nr, np, nsecr))
        wtd <- sweep(ests, MARGIN = 3, STATS = wt, FUN = '*')
        wtd <- apply(wtd, 1:2, sum)
        ## variance from Burnham & Anderson 2004 eq (4)
        vcv1 <- array(unlist(vcvs), dim=c(nr, nr, np, nsecr))  ## between rows of newdata
        betak <- sweep(ests, MARGIN = 1:2, STATS = wtd, FUN='-')
        vcv2 <- array(apply(betak, 2:3, function(x) outer(x,x)), dim=c(nr,nr,np,nsecr))
        vcv <- sweep (vcv1 + vcv2, MARGIN = 4, STATS = wt, FUN = '*')
        vcv <- apply(vcv, 1:3, sum)
        sewtd  <- apply(vcv, 3, function (x) diag(x)^0.5)

    }

    ## OUTPUT
    sewtd <- array(sewtd, dim = c(nr,np))
    z <- abs(qnorm(1-alpha/2))   ## beware confusion with hazard z!
    output <- array (dim=c(nr, 4, np))
    if (type=='beta') {
        # extra dimension just to share later code
        output[1,1,] <- wtd
        output[1,2,] <- sewtd
        output[1,3,] <- wtd - z*sewtd
        output[1,4,] <- wtd + z*sewtd
        ## overwrite confidence intervals if MATA requested
        if (CImethod == 'MATA') {
            for (p in 1:np)
                output[1,3:4,p]  <- MATA(wt, ests[1,p,], estse[1,p,], alpha)
        }
        dimnames(output) <- list(NULL,
            c('beta', 'SE.beta', 'lcl', 'ucl'),
            parnames)
    }
    else { ## type=='real'
        ## 2010 02 15, 2010 02 23, 2015-09-30
        for (m in 1: length(object))
        if (object[[m]]$details$param > 1) {
            stop ("'model.average' cannot handle nonstandard parameterisations")
        }
        for (i in 1:nr) {
            if (average == 'real') {
                output[i,1,] <- wtd[i,]
                output[i,2,] <- sewtd[i,]
                lpwtd <- Xtransform (wtd[i,], links, parnames)
                selpwtd <- se.Xtransform (wtd[i,], sewtd[i,], links, parnames)
                output[i,3,] <- Xuntransform(lpwtd-z*selpwtd, links, parnames)
                output[i,4,] <- Xuntransform(lpwtd+z*selpwtd, links, parnames)
            }
            else {
                output[i,1,] <- Xuntransform(wtd[i,], links, parnames)
                output[i,2,] <- se.Xuntransform(wtd[i,], sewtd[i,], links, parnames)
                output[i,3,] <- Xuntransform(wtd[i,]-z*sewtd[i,], links, parnames)
                output[i,4,] <- Xuntransform(wtd[i,]+z*sewtd[i,], links, parnames)
            }
            ## overwrite confidence intervals if MATA requested
            if (CImethod == 'MATA') {
                for (p in 1:np)
                    output[i,3:4,p] <- untransform(MATA(wt,ests[i,p,], estse[i,p,],alpha), links[[p]])
            }
        }
        dimnames(output) <- list(rownames,
            c('estimate', 'SE.estimate', 'lcl', 'ucl'),
            parnames)
    }

    ## collapse if possible
    if (dim(output)[1] == 1) {
        if (np==1) output <- output[1,,]
        else output <- t(output[1,,])
    }
    else if (np==1) output <- output[,,1]

    if (covar) {
        ## 2011-11-05; 2012-02-16 bug fix
        if (type == 'real') {
            stop("sorry - covar not available at present",
                 " for model-averaged real paramaters")
            dimnames(vcv) <- list (rownames, rownames, parnames)
            output <- list (ma = output, rowvcv = vcv)
        }
        else {
            dimnames(vcv) <- list (parnames, parnames)
            output <- list (ma = output, linkvcv = vcv)
        }
    }
    return(output)
}
############################################################################################

collate <- function (..., realnames = NULL, betanames = NULL, newdata = NULL,
    alpha = 0.05, perm = 1:4, fields = 1:4) {


    ## 2013-04-20
    object <- secrlist(...)
    if (!is.null(names(object)))
        modelnames <- names(object)
    else
        modelnames <- as.character(match.call(expand.dots=FALSE)$...)

    if ( any (!sapply(object, function (x) inherits(x, 'secr'))) )
        stop ("require fitted secr objects")
    if ( length(object) < 2 )
        warning ("only one model")
    if (!is.list(object) | !inherits(object[[1]], 'secr'))
        stop("object must be secr or list of secr")

    type <- 'real'                     ## default

    parnames <- unique(as.vector(unlist(sapply(object,
        function(x) x$realnames))))  ## default

    if (!is.null(realnames))
        parnames <- realnames
    else if (!is.null(betanames)) {
        type <- 'beta'
        parnames <- betanames
    }

    np <- length(parnames)
    nsecr <- length(object)

    ## rudimentary checks for compatible models
    if (nsecr > 1) {
        objnames <- function(i) switch (type,
            real = object[[i]]$realnames, beta = object[[i]]$betanames)
        test <- sapply (2:nsecr, function(i)
            sum(match(parnames, objnames(i), nomatch=0)>0) == np)
        if (!all(test))
            stop ("parameters not found in all models, or incompatible models")
    }

    getLP <- function (object1) {  ## for predicted values of real parameters
        getfield <- function (x) {
            secr.lpredictor (formula = object1$model[[x]], newdata = newdata,
                indx = object1$parindx[[x]], beta = object1$fit$par,
                beta.vcv = object1$beta.vcv, field = x,
                smoothsetup = object1$smoothsetup[[x]])
        }
        sapply (names(object1$model), getfield, simplify = FALSE)
    }

    if (is.null(newdata)) {

        ## form unified 'newdata' containing all necessary predictors

        ## start with list of model-specific newdata --
        ## each component of tempnewdata is a data.frame of newdata
        ## for the corresponding model
        tempnewdata <- lapply (object, secr.make.newdata)
        ## modified 2010 02 14
        column.list <- list(0)
        for (i in 1:nsecr) column.list[[i]] <- as.list(tempnewdata[[i]])
        column.list <- unlist(column.list, recursive = FALSE)

        column.list <- column.list[unique(names(column.list))]
        column.list <- lapply(column.list, unique)
        common <- names(column.list)[names(column.list) %in% names(newdata)]
        column.list[common] <- newdata[common]   ## user input

        ## modified 2010 02 23 to match levels of session covariates to session
        ## not tested with different session covariates in diff sessions
        sessioncovs <- lapply(object, function(x)
            if(!is.null(x$sessioncov)) data.frame(session=session(x$capthist), x$sessioncov)
            else NULL)
        sessioncovs <- sessioncovs[!sapply(sessioncovs, is.null)]
        scn <- as.vector(sapply(sessioncovs, names))
        scn <- match(unique(scn),scn)
        sessioncovs <- as.data.frame(sessioncovs)[,scn]
        sessioncovnames <- unlist(lapply(object, function(x) names(x$sessioncov)))
        sessioncovariate <- names(column.list) %in% sessioncovnames
        newdata <- expand.grid (column.list[!sessioncovariate])
        ## replaced 2010 03 09 to retain sort order
        #   if (nrow(sessioncovs)>0) newdata <- merge(newdata, sessioncovs)
        if (nrow(sessioncovs)>0) {
            for (i in names(sessioncovs)) {
                if (i != 'session') newdata[,i] <- sessioncovs[newdata$session,i]
            }
        }
    }

    nr <- nrow(newdata)
    rownames <- apply(newdata, 1, function(x) paste(names(newdata), '=', x, sep='', collapse=','))
    z <- abs(qnorm(1-alpha/2))   ## beware confusion with hazard z!

    if (type == 'real') {
        predict <- lapply (object, getLP)
        stripped <- lapply(predict, function(x) lapply(x[parnames], function(y) y[, c('estimate','se')] ))
        stripped <- array(unlist(stripped), dim=c(nr, 2, np, nsecr))
    }
    else {
        coefs <- lapply (object, coef)
        stripped <- lapply(coefs, function(x) x[parnames, c('beta','SE.beta')] )
        stripped <- array(unlist(stripped), dim=c(nr, np, 2, nsecr))
        stripped <- aperm(stripped, c(1,3,2,4))
    }
    output <- array (dim=c(nr, 4, np, nsecr))
    if (type=='real') {
        output[,1:2,,] <- stripped
        for (i in 1:nr) {
            for (m in 1:nsecr) {
                ## changed object[[1]]$link to object[[m]]$link following lines 2010 02 14 - seems right
                output[i,1,,m] <- Xuntransform(stripped[i,1,,m, drop=FALSE], object[[m]]$link, parnames)
                output[i,2,,m] <- se.Xuntransform(stripped[i,1,,m, drop=FALSE], stripped[i,2,,m, drop=FALSE], object[[m]]$link, parnames)
                output[i,3,,m] <- Xuntransform(stripped[i,1,,m, drop=FALSE]-z*stripped[i,2,,m, drop=FALSE], object[[m]]$link, parnames)
                output[i,4,,m] <- Xuntransform(stripped[i,1,,m, drop=FALSE]+z*stripped[i,2,,m, drop=FALSE], object[[m]]$link, parnames)
            }
        }
    }
    else { ## type=='beta'
        output[1,1:2,,] <- stripped
        output[1,3,,] <- stripped[1,1,,,drop=FALSE]-z*stripped[1,2,,,drop=FALSE]
        output[1,4,,] <- stripped[1,1,,,drop=FALSE]+z*stripped[1,2,,,drop=FALSE]
    }

    if (type=='real') {
        dimnames(output) <- list(rownames,
            c('estimate', 'SE.estimate', 'lcl', 'ucl'),
            parnames,
            modelnames)
    }
    else {
        dimnames(output) <- list(NULL,
            c('beta', 'SE.beta', 'lcl', 'ucl'),
            parnames,
            modelnames)
    }
    ## default dimensions:
    ## row, model, statistic, parameter
    output <- aperm(output, c(1,4,2,3))

    return(aperm(output[,,fields,,drop=FALSE], perm))

}
############################################################################################

