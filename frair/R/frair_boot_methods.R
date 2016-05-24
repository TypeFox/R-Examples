## frboot methods
print.frboot <- function(x, ...){
    nbootdone <- x$n_boot-x$n_failed
    percsuc <- round(nbootdone/x$n_boot*100,2)
    
    cat('BOOTSTRAPPED FUNCTIONAL RESPONSE FIT\n')
    cat(paste0('\nResponse:            ', x$response))
    cat(paste0('\nDescription:         ', as.character(frair_responses(show=FALSE)[[x$response]][2])))
    cat(paste0('\nOptimised variables: ', paste0(x$optimvars, collapse=', ')))
    cat(paste0('\nFixed variables:     ', ifelse(test=!is.null(x$fixedvars), yes=paste0(x$fixedvars, collapse=', '), no='NA')))
    cat(paste0('\nBootstrap type:      ', ifelse(test=x$stratified, yes='Stratified', no='Ordinary')))
    cat(paste0('\nFit success:         ', percsuc, '% (', nbootdone, ' of ', x$n_boot, ')', sep=''))
    cat(paste0('\nDuplicated fits:     ', x$n_duplicated))
    cat('\n')
    cat('\nCoefficients (original data):\n')
    print(round(x$coefficients,3))
    cat('\n95% BCa confidence intervals (for more info, see ?confint.frboot):\n')
    btconf <- confint(x, citypes='bca')
    
    print(btconf)
    cat('\nNOTE: It is recomended you inspect the raw fit too (see: ?frair_boot)\n')
}

plot.frboot <- function(x, xlab=x$xvar, ylab=x$yvar, ...){
    plot(x$x, x$y, xlab=xlab, ylab=ylab, ...)
}

lines.frboot <- function(x, all_lines=FALSE, tozero=FALSE, bootcol=1, bootalpha=1/sqrt(x$n_boot), ...){
    fitfun <- get(x$response, pos = "package:frair")
    if(tozero){
        zero_answer <- fitfun(0, as.list(x$coefficients))
        if(is.na(zero_answer)){
            warning(c("The supplied function is undefined at zero.\n",
                      "   Plotting to a minimum of 1e-04 instead."))
            lowval <- 1e-04
        } else {
            lowval <- 0
        }
        newx <- seq(from=lowval, to=max(x$x), length.out = 50)
    } else {
        newx <- seq(from=min(x$x), to=max(x$x), length.out = 50)
    }
    
    if(!all_lines){
        # Plot the mean (original) fit
        newy <- fitfun(newx, as.list(x$coefficients))
        lines(newx, newy, ...)
    } else {
        # Plotting bootlines
        bootcol <- adjustcolor(bootcol, alpha.f = bootalpha) # Sort out colour
        bootcoefs <- na.omit(x$bootcoefs)
        outdd <- matrix(ncol=length(newx), nrow=nrow(bootcoefs))
        # Draw the lines
        for(a in 1:nrow(bootcoefs)){
            outdd[a,] <- fitfun(newx, as.list(as.list(bootcoefs[a,])))
        }
        for(a in 1:nrow(outdd)){
            lines(x=newx, y=outdd[a,], col=bootcol, ...)
        }
    }
}

drawpoly <- function(x, upper, lower, ...) UseMethod("drawpoly")

drawpoly.default <- function(x, upper, lower, ...){
    polygon(x=c(x, rev(x), x[1]), y=c(upper, rev(lower), upper[1]), ...)
    # TODO: https://github.com/dpritchard/frair/issues/26
}

drawpoly.frboot <- function(x, ..., probs=c(0.025, 0.975), tozero=FALSE){
    fitfun <- get(x$response, pos = "package:frair")
    if(tozero){
        zero_answer <- fitfun(0, as.list(x$coefficients))
        if(is.na(zero_answer)){
            warning(c("The supplied function is undefined at zero.\n",
                      "   Plotting to a minimum of 1e-04 instead."))
            lowval <- 1e-04
        } else {
            lowval <- 0
        }
        newx <- seq(from=lowval, to=max(x$x), length.out = 50)
    } else {
        newx <- seq(from=min(x$x), to=max(x$x), length.out = 50)
    }
    bootcoefs <- na.omit(x$bootcoefs)
    outdd <- matrix(ncol=length(newx), nrow=nrow(bootcoefs))
    
    cat('\nCalculating polygons.\n\n')
    flush.console()
    for(a in 1:nrow(bootcoefs)){
        outdd[a,] <- fitfun(newx, as.list(as.list(bootcoefs[a,])))
    }
    
    dd <- apply(outdd, 2, quantile, na.rm=T, probs=probs)
    drawpoly(x=newx, upper=dd[2,], lower=dd[1,], ...)
    #polygon(x=c(newx, rev(newx), newx[1]), y=c(dd[1,], rev(dd[2,]), dd[1,1]), ...)
}

confint.frboot <- function(object, parm='all', level=0.95, ..., citypes='all'){
    # Check that the paramters are valid
    if (parm=='all'){
        parm <- object$optimvars
    } else {
        nomatchy <- match(parm, object$optimvars, nomatch = 0) < 1
        if(any(nomatchy)){
            stop('Supplied parameters do not exist!')
        }
    }
    # Check that the types of confidence intervals are valid
    citypesall <- c('norm', 'basic', 'stud', 'perc', 'bca')
    if (citypes=='all'){
        citypes <- citypesall
    } else {
        nomatchy <- match(citypes, citypesall, nomatch = 0) < 1
        if(any(nomatchy)){
            stop('Supplied confidence interval types are not valid!')
        }
    }
    
    runlist <- expand.grid(parm, citypes, stringsAsFactors=FALSE)
    
    # Setup output
    outcis <- list()
    for(a in 1:nrow(runlist)){
        coefname <- runlist[a,1]
        loc <- which(names(object$fit$t0)==coefname)
        locvar <- which(names(object$fit$t0)==paste0(coefname, 'var'))
        type <- runlist[a,2]
        outcis[[coefname]][[type]] <- list()
        bootciout <- fr_catchlist(boot::boot.ci(object$fit, index=c(loc, locvar), conf=level, type=type))
        # A warning or an error
        if(!is.null(bootciout$error)){
            outcis[[coefname]][[type]][['lower']] <- NA
            outcis[[coefname]][[type]][['upper']] <- NA
            outcis[[coefname]][[type]][['bootciout']] <- NA
            outcis[[coefname]][[type]][['errors']] <- bootciout$error
            outcis[[coefname]][[type]][['warnings']] <- NULL
            outcis[[coefname]][[type]][['notes']] <- NULL
        } else {
            # Probably OK
            if(type=='norm'){
                outcis[[coefname]][[type]][['lower']] <- bootciout$val[[4L]][2L]
                outcis[[coefname]][[type]][['upper']] <- bootciout$val[[4L]][3L]
                outcis[[coefname]][[type]][['bootciout']] <- bootciout$val
                outcis[[coefname]][[type]][['errors']] <- NULL
                outcis[[coefname]][[type]][['warnings']] <- bootciout$warnings
                # NB: Normal intervals dont get extra checks
                # Consequently they do not get 'notes'
                outcis[[coefname]][[type]][['notes']] <- NULL
            } else {
                outcis[[coefname]][[type]][['lower']] <- bootciout$val[[4L]][4L]
                outcis[[coefname]][[type]][['upper']] <- bootciout$val[[4L]][5L]
                outcis[[coefname]][[type]][['bootciout']] <- bootciout$val
                outcis[[coefname]][[type]][['errors']] <- NULL
                # Now, some addiitonal checks, a-la print.bootci()
                rng <- range(bootciout$val[[4L]][2:3])
                R <- bootciout$val$R
                # Save this along with other warnings... See ?print.bootci for an explanation
                if((rng[1L] <= 1) || (rng[2L] >= R)){
                    warns <- c(bootciout$warnings,'Warning: Extreme quantiles used. Intervals will be unstable!')
                    outcis[[coefname]][[type]][['warnings']] <- warns
                    outcis[[coefname]][[type]][['notes']] <- NULL
                } else if((rng[1L] <= 10) || (rng[2L] >= R-9)){
                    outcis[[coefname]][[type]][['warnings']] <- bootciout$warnings
                    outcis[[coefname]][[type]][['notes']] <- 'High quantiles used. Intervals may be unstable!'
                }
            }
        }
    }
    class(outcis) <- c('frconf', 'list')
    return(outcis)
}

print.frconf <- function(x, ...){
    # First extract unique errors / warnings / notes
    allinfo <- NULL
    runlist <- NULL
    for(a in 1:length(x)){ # a indexes 'parm'
        for(b in 1:length(x[[a]])){ # b indexes 'citypes'
            errs <- x[[a]][[b]]$errors
            warn <- x[[a]][[b]]$warnings
            note <- x[[a]][[b]]$notes
            if(!is.null(errs)) {allinfo <- rbind(allinfo, data.frame(info=errs, type=1))}
            if(!is.null(warn)) {allinfo <- rbind(allinfo, data.frame(info=warn, type=2))}
            if(!is.null(note)) {allinfo <- rbind(allinfo, data.frame(info=note, type=3))}
            runlist <- rbind(runlist, data.frame(parm=names(x[a]), citype=names(x[[a]][b]), stringsAsFactors=FALSE))
        }
    }
    if(is.null(allinfo)){
        hasnotes <- FALSE
    } else {
        allinfo <- unique(allinfo)
        allinfo <- allinfo[order(allinfo$type),]
        hasnotes <- TRUE
    }
    cinamemap <- c(norm='Normal approx.', basic='Basic', stud='Studentised', perc='Percentile', bca='BCa')
    cat(format('Coefficient', width=13), format('CI Type', width=15), format('Lower', width=8), format('Upper', width=8), if(hasnotes){format('Notes', width=8)}, '\n', sep='')
    for(a in 1:nrow(runlist)){
        dat <- x[[runlist[a,1]]][[runlist[a,2]]]
        info <- c(dat$errors,dat$warnings,dat$notes)
        infoind <- which(allinfo$info %in% info)
        if(length(infoind)==0) {infoind <- '-'}
        
        # Now start printing!
        cat(format(runlist[a,1], width=13), format(cinamemap[runlist[a,2]], width=15), sep='')
        cat(format(ifelse(is.na(dat$lower), '-', as.character(round(dat$lower,3))), width=8), sep='')
        cat(format(ifelse(is.na(dat$upper), '-', as.character(round(dat$upper,3))), width=8), sep='')
        cat(if(hasnotes){paste(infoind, collapse=',')}, '\n', sep='')
        #cat(paste(infoind, collapse=','), '\n', sep='')
        
    }
    if(hasnotes){
        for(a in 1:nrow(allinfo)){
            cat('Note ', a, ': ', as.character(allinfo$info[a]), '\n', sep='')
        }
    }
}

