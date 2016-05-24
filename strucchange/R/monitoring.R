mefp <- function(obj, ...) UseMethod("mefp")

mefp.formula <-
    function(formula, type = c("OLS-CUSUM", "OLS-MOSUM", "RE", "ME", "fluctuation"),
             data=list(), h=1, alpha=0.05, functional = c("max", "range"),
             period=10, tolerance=.Machine$double.eps^0.5,
             CritvalTable=NULL, rescale=NULL, border=NULL, ...)
{
    type <- match.arg(type)
    if(type == "fluctuation") type <- "RE"
    functional <- match.arg(functional)

    histrescale <- rescale
    if(type=="RE") histrescale <- TRUE
    if(is.null(rescale)) {
      if(type=="RE") {
        rescale <- FALSE
	histrescale <- TRUE
      }
      else {
        rescale <- TRUE
	histrescale <- TRUE
      }
    }

    val <- efp(formula, type=type, h=h, data=data, rescale=histrescale)
    val <- mefp(val, alpha=alpha, functional=functional, period=period,
                tolerance=tolerance, CritvalTable=CritvalTable,
                rescale=rescale, border=border, ...)
    if(length(data) == 0)
        val$data <- NULL
    else
        val$data <- deparse(substitute(data))
    val$call <- val$initcall <- match.call()
    return(val)
}

mefp.efp <-
    function(obj, alpha=0.05, functional = c("max", "range"),
             period=10, tolerance=.Machine$double.eps^0.5,
             CritvalTable=NULL, rescale=NULL, border=NULL, ...)
{
    functional <- match.arg(functional)
    if(! (obj$type %in% c("OLS-CUSUM", "OLS-MOSUM", "ME", "RE")))
        stop("efp must be of type `OLS-CUSUM', `OLS-MOSUM', `RE' or `ME'")

    if(is.null(as.list(obj$call)$data)){
       data <- NULL
    }
    else{
       data <- as.character(as.list(obj$call)$data)
    }

    if(is.null(rescale) & (obj$type == "RE")) rescale <- FALSE
    if(is.null(rescale)) rescale <- obj$rescale
    if(is.null(rescale)) rescale <- TRUE

    ## Bonferroni correction
    elemsiglevel <- alpha / obj$nreg
    histcoef <- obj$coefficients
    histsize <- obj$nobs
    winsize <- obj$par
    K <- floor(winsize*obj$nobs)
    sigmahat <- obj$sigma
    Q12 <- obj$Q12

    logPlus <- function(x) ifelse(x<=exp(1),1,log(x))

    switch(obj$type,

    "OLS-CUSUM" = {

        mreSize <- function(a){
            2-2*(pnorm(a)-a*dnorm(a))
        }
        mreCritval <- function(a){
            abs(2*(pnorm(a)-a*dnorm(a))+alpha-2)
        }
        critval <- optimize(mreCritval, c(0,10), tol = tolerance)$minimum

        computeEmpProc <- function(X, y)
        {
          as.vector(cumsum((y - X %*% histcoef)/(sigmahat*sqrt(histsize))))
        }
        computeEstims <- NULL
        functional <- "max"
        if(is.null(border)){
            border <- function(k){
                x <- k/histsize
                sqrt(x*(x-1)*(critval^2 + log(x/(x-1))))
            }
        }
    },

    "OLS-MOSUM" = {
        if(is.null(CritvalTable))
            CritvalTable <- get("monitorMECritvalTable")
        dntab <- dimnames(CritvalTable)
        if(!(winsize %in% dntab[[1]]))
            stop(paste("winsize h =",winsize,"not available, we have:",
                       paste(dntab[[1]], collapse=", ")))
        if(!(period %in% dntab[[2]]))
            stop(paste("period",period,"not available, we have:",
                       paste(dntab[[2]], collapse=", ")))
        critval <- approx(x=as.numeric(dntab[[3]]),
                          y=CritvalTable[as.character(winsize),
                          as.character(period),,functional],
                          xout=1-alpha)$y
        if(is.na(critval))
            stop(paste("Necessary significance level per parameter of",
                       alpha,
                       "\n\toutside of available range",
                       paste(range(1-as.numeric(dntab[[3]])),
                             collapse="-")))

        computeEmpProc <- function(X, y)
        {
               e <- as.vector(y - X %*% histcoef)
               process <- rep(0, nrow(X)-K+1)
               for(i in 0:(nrow(X)-K))
               {
                   process[i+1] <- sum(e[(i+1):(i+K)])
               }
               process/(sigmahat*sqrt(histsize))
        }
        functional <- "max"
        computeEstims <- NULL
        if(is.null(border)){
            border <- function(k){
                critval*sqrt(2*logPlus(k/histsize))
            }
        }
    },

    "ME" = {
        if(is.null(CritvalTable))
            CritvalTable <- get("monitorMECritvalTable")
        dntab <- dimnames(CritvalTable)
        if(!(winsize %in% dntab[[1]]))
            stop(paste("winsize h =",winsize,"not available, we have:",
                       paste(dntab[[1]], collapse=", ")))
        if(!(period %in% dntab[[2]]))
            stop(paste("period",period,"not available, we have:",
                       paste(dntab[[2]], collapse=", ")))
        critval <- approx(x=as.numeric(dntab[[3]]),
                          y=CritvalTable[as.character(winsize),
                          as.character(period),,functional],
                          xout=1-elemsiglevel)$y
        if(is.na(critval))
            stop(paste("Necessary significance level per parameter of",
                       elemsiglevel,
                       "\n\toutside of available range",
                       paste(range(1-as.numeric(dntab[[3]])),
                             collapse="-")))

        computeEstims <- function(x, y, k){
            ok <- (k-K+1):k
            retval <- list(coef=NULL, Qr12=NULL)
            retval$coef <- coef(lm.fit(x[ok,,drop=FALSE], y[ok,,drop=FALSE]))
            if(rescale)
                retval$Qr12 <- root.matrix(crossprod(x[ok,,drop=FALSE]))/sqrt(K)
            retval
        }

        computeEmpProc <- function(newcoef, Q, k){
            if(is.null(Q)) Q <- Q12
            Q <- Q * K/(sigmahat*sqrt(histsize))
            t(Q %*%(newcoef-histcoef))
        }

        if(is.null(border)){
            border <- function(k){
                critval*sqrt(2*logPlus(k/histsize))
            }
        }
    },

    "RE" = {

        if(is.null(CritvalTable)){
            mreSize <- function(a){
                2-2*(pnorm(a)-a*dnorm(a))
            }
            mreCritval <- function(a){
                abs(2*(pnorm(a)-a*dnorm(a))+elemsiglevel-2)
            }
            critval <- optimize(mreCritval, c(0,10), tol = tolerance)$minimum
        }
        else{
            dntab <- dimnames(CritvalTable)
            if(!(period %in% dntab[[1]]))
                stop(paste("period",period,"not available, we have:",
                           paste(dntab[[1]], collapse=", ")))
            critval <- approx(x=as.numeric(dntab[[2]]),
                              y=CritvalTable[as.character(period),],
                              xout=1-elemsiglevel)$y
            if(is.na(critval))
                stop(paste("Necessary significance level per parameter of",
                           elemsiglevel,
                           "\n\toutside of available range",
                           paste(range(1-as.numeric(dntab[[3]])),
                                 collapse="-")))
        }


        computeEstims <- function(x, y, k){
            retval <- list(coef=NULL, Qr12=NULL)
            retval$coef <- coef(lm.fit(x[1:k,,drop=FALSE], y[1:k,,drop=FALSE]))
            if(rescale)
                retval$Qr12 <- root.matrix(crossprod(x[1:k,,drop=FALSE]))/sqrt(k)
            retval
        }

        computeEmpProc <- function(newcoef, Q, k){
            if(is.null(Q)) Q <- Q12
            Q <- Q * k/(sigmahat*sqrt(histsize))
            t(Q %*%(newcoef-histcoef))
        }

        if(is.null(border)){
            border <- function(k){
                x <- k/histsize
                sqrt(x*(x-1)*(critval^2 + log(x/(x-1))))
            }
        }
    })

    if(functional=="max"){
        computeStat <- function(empproc){
            max(abs(empproc))
        }
    }
    else if(functional=="range"){
        if(obj$type=="RE")
            stop("Functional `range' not available for recursive estimates")
        else{
            computeStat <- function(empproc){
                max(apply(empproc, 2, function(x) diff(range(x))))
            }
        }
    }

    obj <- list(breakpoint=NA, last=obj$nobs, process=NULL,
                statistic=NULL, histsize=histsize,
                initcall=match.call(), call=match.call(),
                efpcall=obj$call, efpprocess=obj$process,
                computeEstims=computeEstims,
                computeEmpProc=computeEmpProc,
                border=border, computeStat=computeStat,
                functional=functional, alpha=alpha, critval=critval,
                histcoef=histcoef, formula=obj$formula,
                type.name=paste("Monitoring with", obj$type.name),
                type=obj$type, data=data, histtsp=obj$datatsp)

    class(obj) <- "mefp"
    obj
}

monitor <- function(obj, data=NULL, verbose=TRUE){

    if(!is.na(obj$breakpoint)) return(TRUE)
    if(missing(data)){
        if(is.null(obj$data)){
            data <- list()
        }
        else{
            data <- get(obj$data)
        }
    }

    mf <- model.frame(obj$formula, data=data)
    y <- as.matrix(model.response(mf))
    modelterms <- terms(obj$formula, data = data)
    x <- model.matrix(modelterms, data = data)

    if(nrow(x) <= obj$last) return(obj)
    if(nrow(x)!=nrow(y))
        stop("response and regressors must have the same number of rows")
    if(ncol(y)!=1)
        stop("multivariate response not implemented yet")
    foundBreak <- FALSE

    if((obj$type == "OLS-MOSUM") | (obj$type == "OLS-CUSUM"))
    {
      if(obj$type == "OLS-CUSUM")
      {
        obj$process <- obj$computeEmpProc(x,y)[-(1:obj$histsize)]
      }
      else
      {
        obj$process <- obj$computeEmpProc(x,y)[-(1:length(obj$efpprocess))]
      }
      boundary <- obj$border((obj$histsize+1):nrow(x))
      obj$statistic <- max(abs(obj$process))
      if(!foundBreak & any(abs(obj$process) > boundary))
      {
        foundBreak <- TRUE
        obj$breakpoint <- min(which(abs(obj$process) > boundary)) + obj$histsize
        if(verbose) cat("Break detected at observation #", obj$breakpoint, "\n")
      }
      obj$lastcoef <- NULL
    }
    else
    {
      for(k in (obj$last+1):nrow(x)){
          newestims <- obj$computeEstims(x,y,k)
          obj$process <- rbind(obj$process,
                               obj$computeEmpProc(newestims$coef, newestims$Qr12,k))
          stat <- obj$computeStat(obj$process)
          obj$statistic <- c(obj$statistic, stat)
          if(!foundBreak & (stat > obj$border(k))){
              foundBreak <- TRUE
              obj$breakpoint <- k
              if(verbose) cat("Break detected at observation #", k, "\n")
          }
      }
      obj$lastcoef <- newestims$coef
    }
    obj$last <- nrow(x)
    obj$call <- match.call()
    obj
}

print.mefp <- function(x, ...){

    cat(x$type.name, "\n\n")
    cat("Initial call:\n ", deparse(x$initcall), "\n\n")
    cat("Last call:\n ", deparse(x$call), "\n\n")
    cat("Significance level   : ", x$alpha, "\n")
    cat("Critical value       : ", x$critval, "\n")
    cat("History size         : ", x$histsize, "\n")
    cat("Last point evaluated : ", x$last, "\n")
    if(!is.na(x$breakpoint))
        cat("Structural break at  : ", x$breakpoint, "\n")
    cat("\nParameter estimate on history :\n");
    print(x$histcoef)
    if(!is.null(x$lastcoef)){
        cat("Last parameter estimate :\n");
        print(x$lastcoef)
    }
}

plot.mefp <- function(x, boundary=TRUE, functional="max", main=NULL,
                      ylab="Empirical fluctuation process", ylim=NULL, ...){

    if(x$last>x$histsize){
        proc <- rbind(as.matrix(x$efpprocess),
                    as.matrix(x$process))
        proc <- ts(proc,
                   end=x$histtsp[2]+(x$last-x$histsize)/x$histtsp[3],
                   frequency=x$histtsp[3])
        bound <- ts(x$border((x$histsize+1):x$last),
                 end = end(proc), frequency=frequency(proc))
        pos <- FALSE
        if(!is.null(functional) && (functional == "max"))
        {
            proc <- ts(apply(abs(proc), 1, 'max'),
                       start = start(proc), frequency = frequency(proc))
            pos <- TRUE
        }
        ymax <- max(c(proc, bound))
        if(pos)
            ymin <- 0
        else
            ymin <- min(c(proc, -bound))
        if(is.null(ylim)) ylim <- c(ymin, ymax)
        if(is.null(main))
            main <- x$type.name
        if(boundary)
            panel <- function(y, ...)
            {
                lines(y, ...)
                lines(bound, col=2)
                lines(-bound, col=2)
                abline(0,0)
            }
        else
            panel <- function(y, ...)
            {
                lines(y, ...)
                abline(0,0)
            }
        if(any(attr(proc, "class") == "mts"))
            plot(proc, main = main, panel = panel, ...)
        else
        {
            plot(proc, main = main, ylab = ylab, ylim = ylim, ...)
            if(boundary)
            {
                lines(bound, col=2)
                if(!pos) lines(-bound, col=2)
            }
            abline(0,0)
        }
        abline(v=x$histtsp[2], lty=2)
    }
    else{
        cat("Nothing monitored yet!\n")
    }
}

boundary.mefp <- function(x, ...)
{
    ts(x$border((x$histsize+1):x$last),
       start = x$histtsp[2]+1/x$histtsp[3],
       frequency=x$histtsp[3])
}

lines.mefp <- function(x, ...)
{
    if(x$last>x$histsize){
        proc <- rbind(as.matrix(x$efpprocess),
                    as.matrix(x$process))
        proc <- ts(proc,
                   end=x$histtsp[2]+(x$last-x$histsize)/x$histtsp[3],
                   frequency=frequency(x$efpprocess))
        proc <- ts(apply(abs(proc), 1, 'max'),
                   start = start(proc), frequency = frequency(proc))
        lines(proc, ...)
    }
    else{
        cat("Nothing monitored yet!\n")
    }
}


