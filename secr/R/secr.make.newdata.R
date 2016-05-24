############################################################################################
## package 'secr'
## secr.make.newdata.R
## last changed
## 2009 12 13 (mixtures)
## 2010 03 10 'T'
## 2010 06 17 'Session'
## 2010 06 21 'x2', 'y2', 'xy'
## 2010 08 28 fix bug with T
## 2011 11 28 user dframe factors now covered
## 2015-10-08 'ts'
## Create (neutral) design data suitable for 'predict'
############################################################################################

secr.make.newdata <- function (object) {

    findvars.MS <- function (cov, vars, dimcov, use.all) {
        ## function to add covariates to a design data frame 'dframe'
        ## cov may be a dataframe or list of dataframes, one per session (R > 1),
        ## if list, then require predictors to appear in all sessions
        ## uses pad1 and insertdim from functions.R
        ## NOT to be used to add group variables
        ## Does not yet standardize numeric covariates if (!is.factor(vals)) vals <- scale(vals)

        if (is.null(cov) | (length(cov)==0) | (length(vars)==0)) return()
        else {
            found <- ''
            if (!is.data.frame(cov)) {   ## therefore assume is a list
                if (!is.list(cov) | (R==1))
                    stop ('irregular covariates - check multisession structure')
                covnames <- lapply(cov, names)
                varincov <- sapply(covnames, function(nam) vars %in% nam)
                if (length(vars)>1) found <- vars[apply(varincov,1,all)]
                else found <- vars[all(varincov)]

                for (variable in found) {
                    ## use first occurrence!
                    vals <- unlist(lapply(cov, function(x) rep(x[1,variable],
                        dims[dimcov])))
                    newdata[,variable] <<- insertdim (vals, dimcov, dims)
                }
            }
            else
            {
                found <- names(cov) %in% vars
                if (is.data.frame(cov) & any(found)) {
                    found <- names(cov)[found]
                    values <- as.data.frame(cov[,found])
                    names(values) <- found
                    if (length(values)>0) {
                        for (variable in found) {
                            if (use.all) vals <- values[,variable]
                            else  vals <- values[1,variable]
                            newdata[,variable] <<- insertdim (vals, dimcov, dims)
                        }
                    }
                }
            }
            vars <<- vars[!(vars %in% found)]
        }
    }

    capthist <- object$capthist
    mask <- object$mask
    vars <- object$vars
    groups <- object$groups
    timecov <- object$timecov
    sessioncov <- object$sessioncov
    nmix <- object$details$nmix
    hcov <- object$hcov

    if(is.null(nmix)) nmix <- 1
    mixvar <- switch(nmix, character(0),'h2','h3')

    nocc <- max(n.occasion (capthist))
    grouplevels <- group.levels(capthist, groups)
    ngrp <- max(1, length(grouplevels))
    sessions <- session(capthist)
    R <- length(sessions)
    dims <- c(R, ngrp, nmix)

    basevars <- list(session = sessions)
    if (ngrp>1) basevars$g <- factor(grouplevels)
    if (nmix>1) basevars[mixvar] <- list(h.levels(capthist, hcov, nmix))
    newdata <- expand.grid(basevars)
    nr <- nrow(newdata)  ## one row for each session, group and mixture
    if (ngrp==1)
        findvars.MS (covariates(capthist), vars, 1, FALSE) ## check for indiv cov
    for (v in vars) {
        if (v=='x') newdata$x <- rep(0,nr)   # mean attr(mask,'meanSD')[1,'x']
        if (v=='y') newdata$y <- rep(0,nr)   # mean attr(mask,'meanSD')[1,'y']
        if (v=='x2') newdata$x2 <- rep(0,nr)   # mean attr(mask,'meanSD')[1,'x']
        if (v=='y2') newdata$y2 <- rep(0,nr)   # mean attr(mask,'meanSD')[1,'y']
        if (v=='xy') newdata$xy <- rep(0,nr)   # mean attr(mask,'meanSD')[1,'x']
        if (v=='t') newdata$t <- rep(factor(1, levels=1:nocc), nr)   ## mod 2009 09 03
        if (v=='ts') newdata$ts <- rep(factor('marking', levels=c('marking','sighting')), nr)   ## 2015-10-08
        if (v=='T') newdata$T <- rep(0, nr)   ## 2010 08 28
        if (v=='b') newdata$b <- rep(factor(0, levels=c(0,1)),nr)    # naive
        if (v=='B') newdata$B <- rep(factor(0, levels=c(0,1)),nr)    # naive
        if (v=='bk') newdata$bk <- rep(factor(0, levels=c(0,1)),nr)   # naive
        if (v=='Bk') newdata$Bk <- rep(factor(0, levels=c(0,1)),nr)   # naive
        if (v=='k') newdata$k <- rep(factor(0, levels=c(0,1)),nr)    # naive
        if (v=='K') newdata$K <- rep(factor(0, levels=c(0,1)),nr)    # naive
        if (v=='bkc') newdata$bkc <- rep(factor('None', levels=c('None','Self','Other','Both')),nr)
        if (v=='Bkc') newdata$Bkc <- rep(factor('None', levels=c('None','Self','Other','Both')),nr)
        if (v=='tcov') {
            timecov <- object$timecov
            if (is.factor(timecov)) {
                newdata$tcov <- rep(factor(levels(timecov)[1], levels = levels(timecov)))
            }
            else
                newdata$tcov <- rep(0,nr)        # ideally use mean or standardize?
        }
        if (v=='kcov') {
            kcov <- covariates(traps(object$capthist))[,1]
            if (is.factor(kcov)) {
                newdata$kcov <- rep(factor(levels(kcov)[1], levels = levels(kcov)))
            }
            else {
                newdata$kcov <- rep(0,nr)        # ditto
            }
        }
        if (v=='Session') newdata$Session <- as.numeric( factor(newdata$session,
            levels = session(capthist) ) ) - 1    # based on sequence in capthist
    }

    ## all autovars should now have been dealt with
    vars <- vars[!vars %in% c('g','x','y','x2','y2','xy','session','Session',
        't','T','ts','b','B','bk','Bk','bkc','Bkc','k','K','tcov','kcov','h2','h3')]

    findvars.MS (sessioncov, vars, 1, TRUE)
    findvars.MS (timecov, vars, 1, FALSE)
    findvars.MS (covariates(traps(capthist)), vars, 1, FALSE)
    ## added 2011-11-14
    findvars.MS (covariates(mask), vars, 1, FALSE)

    ## 2011-11-28
    if (!is.null(object$dframe)) {
        dframevars <- names(object$dframe)
        for (v in dframevars)
            newdata[,v] <- rep(object$dframe[1,v], nr)
        vars <- vars[!vars %in% dframevars]
    }

    ## default all remaining vars to numeric zero
    for (v in vars) newdata[,v] <- rep(0,nr)

    newdata
}
############################################################################################

