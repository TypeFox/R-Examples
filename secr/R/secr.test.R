############################################################################################
## package 'secr'
## secr.test.R
## last changed 2014-08-07, 2014-09-07
############################################################################################

## Monte Carlo goodness of fit tests
## based either on a statistic computed on raw data without re-fitting model (fit = FALSE)
## or on a statistic computed on each fitted model (fit = TRUE)
## default: proportion detected more than once, fit = FALSE

## query: does this simulate single-catch if object$capthist is single?
## need to make safe for ms
## if ( detector(traps(object$capthist)) == 'single')
##    detector(traps(object$capthist)) <- 'multi'

############################################################################################

secr.test <- function (object, nsim = 99, statfn, fit = FALSE, seed = NULL,
                       ncores = 1, tracelevel = 1) {

    summarise <- function (CH.or.secr, sims1) {
        observed <- statfn(CH.or.secr)
        simulated <- sapply(sims1, statfn)
        nstat <- length(observed)
        rown <- names(observed)
        if (is.null(rown)) {
            rown <- paste('statistic', 1:nstat, sep = '')
            names(observed) <- rown
        }
        observedmat <- matrix(observed, nrow = nstat, ncol = 1, dimnames = list(rown,NULL))
        simulated <- matrix(simulated, nrow = nstat, ncol = nsim, dimnames = list(rown,NULL))
        ptot <- cbind( observedmat, simulated)
        p <- apply(ptot, 1, function(x) rank(x)[1]/ (nsim+1))
        list(simulated = simulated, observed = observed, p = p)
    }

    sims <- simulate(object = object, nsim = nsim, ncores = ncores, seed = seed)

    if (fit) {
        if (missing(statfn))
            ## object is secr fit
            statfn <- function(object) c(devdf = deviance(object) / df.residual(object))
        fitted <- sim.secr(object, nsim = nsim, extractfn = trim, seed = seed,
                           data = sims, start = object$fit$par, ncores = ncores,
                           tracelevel = tracelevel)
        out <- summarise(object, fitted)
    }
    else {
        if (missing(statfn))
            ## object is capthist
            statfn <- function(object) c(f1 = sum(apply(abs(object) > 0, 1, sum) == 1) /
                                         nrow(object))
        if (ms(object)) {
            ## turn sims inside out
            nsess <- length(object$capthist)
            newsims <- vector(nsess, mode = 'list')
            ## reorganise simulations into list by session, each new
            ## component is a list of nsims single-session capthist
            for (j in 1:nsess)
                newsims[[j]] <- lapply(sims, '[[', j)
            out <- mapply(summarise, object$capthist, newsims, SIMPLIFY = FALSE)
        }
        else {
            out <- summarise(object$capthist, sims)
        }
    }

    output <- list(object = object, nsim = nsim, statfn = statfn,
                   fit = fit, seed = seed, output = out)
    class(output) <- c('secrtest', 'list')
    output
}
############################################################################################

ms.secrtest <- function(object, ...) {
    !object$fit & ms(object$object)
}
############################################################################################

plot.secrtest <- function(x, stat, ...) {
    plotcomp <- function(y) {
        nstat <- length(stat)
        rown <- names(y$observed)
        if (!all(stat %in% rown))
            stop ("error in stat name(s)")
        plotone <- function(i) {
            out <- hist(y$simulated[i,], xlab = i, ...)
            abline (v = y$observed[i], col = 'red')
            out
        }
        sapply(stat, plotone)
    }
    if (ms(x)) {
        if (missing(stat))
            stat <- names(x$output[[1]]$observed)
        invisible(lapply(x$output, plotcomp))
    }
    else {
        if (missing(stat))
            stat <- names(x$output$observed)
        invisible(plotcomp(x$output))
    }
}
############################################################################################

print.secrtest <- function (x, terse = TRUE, ...) {
    if (ms(x)) {
        if (terse) {
            getp <- lapply(x$output, '[[', 'p')
            nsess <- length(getp)
            getp <- matrix (unlist(getp), ncol = nsess, byrow = TRUE,
                            dimnames= list(names(x$output[[1]]$observed),
                                    names(getp)))
            print(getp)
        }
        else {
            print(c(x[-6],x['output']))
        }
    }
    else {
        if (terse) {

            print(data.frame(p = x$output$p))
        }
        else {
            NextMethod('print', x)  ## print.list
        }
    }
}
############################################################################################

## Freeman-Tukey statistic?
##

## Why is simulate.secr so slow?
## > system.time(simulate(secrdemo.0, nsim=100))
##    user  system elapsed
##   52.09    0.71   52.90
## > system.time(for(i in 1:100) sim.capthist(traps(captdata), detectpar = detectpar(secrdemo.0)))
##    user  system elapsed
##    1.49    0.00    1.48

##  Rprof()
## > tmp <- simulate(secrdemo.0, nsim=100)
## > Rprof(NULL)
## > summaryRprof()
##
## ...
## $by.total
##                         total.time total.pct self.time self.pct
## "simulate"                   50.42    100.00      0.00     0.00
## "simulate.secr"              50.42    100.00      0.00     0.00
## "FUN"                        50.40     99.96      0.22     0.44
## "lapply"                     50.40     99.96      0.02     0.04
## "sim.detect"                 50.26     99.68      0.56     1.11
## "secr.design.MS"             47.62     94.45      0.04     0.08
## "sapply"                     34.34     68.11      0.02     0.04
## ...

## so secr.design.MS is the problem...
