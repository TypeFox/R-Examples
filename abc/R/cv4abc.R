######################################################################
#
# cv4abc.R
#
# copyright (c) 2011-05-30, Katalin Csillery, Olivier Francois and
# Michael GB Blum
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
# 
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
# 
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
# 
# Part of the R/abc package
# Contains: cv4abc, is.cv4abc, summary.cv4abc, plot.cv4abc
#
######################################################################


cv4abc <- function(param, sumstat, abc.out = NULL, nval, tols,
                   statistic = "median", prior.range = NULL,
                   method, hcorr = TRUE,
                   transf = "none", logit.bounds = c(0,0), subset = NULL, kernel = "epanechnikov",
                   numnet = 10, sizenet = 5, lambda = c(0.0001,0.001,0.01), trace = FALSE, maxit = 500, ...){

    mywarn <- options()$warn
    options(warn=-1)
    linout <- TRUE
    ## checks:
    ## ########
    if(!any(statistic == c("median", "mean", "mode"))){
        stop("Statistic has to be mean, median or mode.", call.=F)
    }
    if(is.null(abc.out) && missing(method)) stop("Method must be supplied when 'abc.out' is NULL.", call.=F)
    if(missing(nval)) stop("'nval' must be supplied.", call.=F)
    
    ## set random seeds
    ## ################
    if(!exists(".Random.seed", envir=.GlobalEnv, inherits = FALSE)) runif(1)
    seed <- get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)
    
    ## define defaults:
    ## #################
    
    if(!is.null(abc.out)){
        subset <- abc.out$na.action
        method <- abc.out$method
        transf <- abc.out$transf
        logit.bounds <- abc.out$logit.bounds
        kernel <- "epanechnikov"
    }
    
    ## checks, numbers of stats, params, sims
    if(is.null(dim(param))){
        np <- 1
        param <- as.data.frame(param)
        names(param) <- "P1"
    }
    else np <- dim(param)[2]
    if(is.null(dim(sumstat))){
        numstat <- 1
        sumstat <- as.data.frame(sumstat)
        names(sumstat) <- "S1"
    }
    else{
        numstat <- dim(sumstat)[2]
    }
    numsim <- dim(sumstat)[1]
    
    ## paramnames & statnames 
    if(!is.null(abc.out)){ # paramnames & statnames from abc.out
        if(np != abc.out$numparam || numstat != abc.out$numstat || numsim != length(abc.out$na.action)){
            stop("The number of parameters, summary statistics, or simulations provided in 'param' or 'sumstat' are not the same as in 'abc.out'.", call.=F)
        }
        else if(!prod(colnames(param) %in% abc.out$names$parameter.names)){
            stop("Parameters in 'param' are not the same as in 'abc.out', or different names are used.", call.=F)
        }
        else if(!prod(colnames(sumstat) %in% abc.out$names$statistics.names)){
            stop("Summary statistics in 'sumstat' are not the same as in 'abc.out', or different names are used.", call.=F)
        }
        else{
            paramnames <- abc.out$names$parameter.names
            statnames <- abc.out$names$statistics.names
        }
    }
    else{ # give paramnames & statnames o/w
        if(length(colnames(param))){
            paramnames <- colnames(param)
        }
        else{
            paramnames <- paste("P", 1:np, sep="")
        }
    }
    if(length(colnames(sumstat))){
        statnames <- colnames(sumstat)
    }
    else{
        statnames <- paste("S", 1:numstat, sep="")
    }
    
    ## indexes for the CV sample and check that the sample is not actually an NA
    gwt <- rep(TRUE,length(sumstat[,1]))
    gwt[attributes(na.omit(sumstat))$na.action] <- FALSE
    if(is.null(subset)) subset <- rep(TRUE,length(sumstat[,1]))
    gwt <- as.logical(gwt*subset)
    cvsamp <- sample(1:numsim, nval, prob = gwt/sum(gwt))
    
    ## if tolerances is a vector have to loop through all values
    
    ## Do all this stuff one by one for the parameteres
    tols <- sort(tols)
    num.panel <- length(tols)
    alltol <- list()
    mycall <- list()
    for(mytol in tols){
        res <- matrix(ncol = np, nrow = nval)
        for(i in 1:nval){
            ## things to over-write from original call: tolerances, target, param, sumstat
            mysamp <- cvsamp[i]
            mytrue <- param[mysamp,]
            mytarget <- sumstat[mysamp,]
            myparam <- param[-mysamp,]
            mysumstat <- sumstat[-mysamp,]
            mysubset <- subset[-mysamp]
            subres <- withCallingHandlers( abc(target = mytarget, param = myparam, sumstat = mysumstat, tol=mytol,
                          subset = mysubset, method = method, transf = transf, logit.bounds = logit.bounds, kernel = kernel,hcorr = hcorr), warning = namesWarningFilter)
            if(statistic == "median") estim <- invisible(summary.abc(subres, print = F, ...)[3, ])
            if(statistic == "mean") estim <- invisible(summary.abc(subres, print = F, ...)[4, ])
            if(statistic == "mode") estim <- invisible(summary.abc(subres, print = F, ...)[5, ])
            res[i, ] <- estim
        }
        if(np == 1) res <- c(res)
        else colnames(res) <- paramnames
        alltol[[paste("tol", mytol, sep="")]] <- res
        mycall[[paste("tol", mytol, sep="")]] <- call("abc", target = quote(target), param = quote(param), sumstat = quote(sumstat),
                                                      tol= mytol,
                                                      subset = quote(subset), method = subres$method, hcorr = subres$hcorr,transf = subres$transf,
                                                      logit.bounds = subres$logit.bounds, kernel = subres$kernel)
    }
    if(np==1){
        true <- as.data.frame(param[cvsamp,])
        names(true) <- paramnames
    }
    else true <- param[cvsamp,]
    cv4abc.out <- list(calls = mycall, cvsamples = cvsamp, tols = tols, true = true, estim = alltol,
                       names = list(parameter.names=paramnames, statistics.names=statnames), seed=seed)

    options(warn=mywarn)
    class(cv4abc.out) <- "cv4abc"
    invisible(cv4abc.out)
    
}

is.cv4abc <- function(x){
    if (inherits(x, "cv4abc")) TRUE
    else FALSE    
}

plot.cv4abc <- function(x, exclude = NULL, log = NULL, file = NULL, postscript = FALSE, onefile = TRUE, ask = !is.null(deviceIsInteractive()), caption = NULL, ...){
    
    if (!inherits(x, "cv4abc")) 
        stop("Use only with objects of class \"cv4abc\".", call.=F)
  
    cv4abc.out <- x
    tols <- cv4abc.out$tols
    numtols <- length(tols)
    np <- length(cv4abc.out$names$parameter.names)
    true <- cv4abc.out$true
    estim <- cv4abc.out$estim
    parnames <- cv4abc.out$names$parameter.names
    cv4abc.out$estim <- as.data.frame(cv4abc.out$estim)
    nval <- length(true)/np
    
    if(is.null(log)) log <- rep("", each=np)
    else if(length(log) != np) stop("error in argument 'log': provide scale for all parameters.", call.=F)
    
    if(is.null(caption)) caption <- as.graphicsAnnot(parnames)
    
    cols <- heat.colors(numtols)
    pch <- 20
    
    ## Devices
    ## ##########
    save.devAskNewPage <- devAskNewPage()
    if(!is.null(file)){
        file <- substitute(file)
        if(!postscript) pdf(file = paste(file, "pdf", sep="."), onefile=onefile)
        if(postscript) postscript(file = paste(file, "ps", sep="."), onefile=onefile)
    }
    else{
        if (ask && 1 < np) {
            devAskNewPage(TRUE)
        }
    }
    
    par(cex = 1, cex.main = 1.2, cex.lab = 1.1)
    for(par in 1:np){
        
        mylog <- log[par]
        columns <- seq(par, numtols*np, by=np)
        
        if(!is.null(exclude)) plot(rep(cv4abc.out$true[-exclude,par], times=numtols),
                                   unlist(cv4abc.out$estim[-exclude,columns]),
                                   col = rep(cols, each=nval),
                                   pch = pch, log=mylog,
                                   xlab="True value", ylab="Estimated value", main=caption[par])
        else  plot(rep(cv4abc.out$true[,par], times=numtols),
                   unlist(cv4abc.out$estim[,columns]),
                   col = rep(cols, each=nval),
                   pch = pch, log=mylog,
                   xlab="True value", ylab="Estimated value", main=caption[par])
        
        abline(0,1)
    }
    
    if(!is.null(file)){
        dev.off()
    }
    else devAskNewPage(save.devAskNewPage)
    
    invisible()
    
}

summary.cv4abc <- function(object, print = TRUE, digits = max(3, getOption("digits")-3), ...){
    
    if (!inherits(object, "cv4abc")) 
        stop("Use only with objects of class \"cv4abc\".", call.=F)
    
    cv4abc.out <- object
    tols <- cv4abc.out$tols
    numtols <- length(tols)
    np <- length(cv4abc.out$names$parameter.names)
    true <- cv4abc.out$true
    estim <- cv4abc.out$estim
    
    parnames <- cv4abc.out$names$parameter.names
    nval <- length(cv4abc.out$cvsamples)
    
    cat(paste("Prediction error based on a cross-validation sample of ",nval,"\n\n", sep=""))
    
    ## error for all parameters and tolerances
    sqdiff <- lapply(estim, function(a) apply((a-true)^2, 2, sum))
    truevar <- apply(true, 2, var)*nval
    prederr <- lapply(sqdiff, function(a) a/truevar)
    prederr <- t(as.data.frame(prederr))
    rownames(prederr) <- tols
    
    class(prederr) <- "table"
    prederr
    
}
