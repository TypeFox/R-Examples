###############################################################################
## package 'secrdesign'
## methods.R
## 2014-02-07,08
## 2014-04-09, 2014--04-11, 2014-04-29
## 2014-11-25 select.stats moved to own file
## 2015-11-26 collapsesubarg for header()
###############################################################################

make.array <- function (object) {
    inputs <- attr(object$scenarios, 'inputs')
    if (is.null(inputs))
        stop("array output requires 'inputs' attribute")
    dims <- sapply(inputs, length)
    vdims <- dims[dims>1]
    varying <- names(inputs)[dims > 1]
    statistics <- dimnames(object$output[[1]])[[2]]
    nstat <- length(statistics)
    nrepl <- nrow(object$output[[1]])
    if (length(varying)>0) {
        outputarray <- array (dim = c(nrepl, nstat, vdims), dimnames =
                              c(list(1:nrepl), list(statistics), inputs[varying]))
    }
    else {
        outputarray <- array (dim = c(nrepl, nstat), dimnames =
                              list(NULL, statistics))
    }
    outputarray[] <- unlist(object$output)
    outputarray
}

predict.fittedmodels <- function (object, ...) {
    output <- lapply(object$output, lapply, predict, ...)
    object$output <- output
    object$outputtype <- 'predicted'
    class(object) <- c('estimatetables', 'secrdesign', 'list')
    object
}

coef.fittedmodels <- function (object, ...) {
    output <- lapply(object$output, lapply, coef, ...)
    object$output <- output
    object$outputtype <- 'coef'
    class(object) <- c('estimatetables', 'secrdesign', 'list')
    object
}

derived.SL <- function (object, ...) {
    if (!inherits(object,'fittedmodels'))
        stop ("require fitted secr models")
    output <- lapply(object$output, lapply, derived, ...)
    object$output <- output
    object$outputtype <- 'derived'
    class(object) <- c('estimatetables', 'secrdesign', 'list')
    object
}

regionN.SL <- function (object, ...) {
    if (!inherits(object,'fittedmodels'))
        stop ("require fitted secr models")
    output <- lapply(object$output, lapply, region.N, ...)
    ra <- sapply(output, function(x) attr(x[[1]], 'regionsize'))
    object$output <- output
    object$outputtype <- 'regionN'
    attr(object, 'regionsize') <- ra    ## one per scenario
    class(object) <- c('estimatetables', 'secrdesign', 'list')
    object
}

## need to handle args that are function or large userdist matrix
## 2014-11-12, 23

dummyuserdist <- function (arg) {
    if (is.list(arg)) {
        if (!is.null(arg$details$userdist)) {
            if (is.function(arg$details$userdist))
                arg$details$userdist <- 'userdistfn'
            else
                arg$details$userdist <- 'userdistmatrix'
        }
        if (!is.null(arg$userdist)) {
            if (is.function(arg$userdist))
                arg$userdist <- 'userdistfn'
            else
                arg$userdist <- 'userdistmatrix'
        }
    }
    arg
}

collapsemodel <- function (arg) {
    if ('model' %in% names(arg)) {
        arg[['model']] <-  secr:::model.string(arg[['model']], NULL)
    }
    arg
}

## 2015-11-26
collapsesubarg <- function (arg,subarg) {
    if (subarg %in% names(arg)) {
        ab <- arg[[subarg]]
        arg[[subarg]] <- paste(paste(names(ab),ab,sep='='), collapse=',')
    }
    arg
}

argdf <- function (args) {
    if (is.null(args))
        NULL
    else { 
        tmp <- lapply(args, dummyuserdist)
        tmp <- lapply(tmp, collapsemodel)
        tmp <- lapply(tmp, collapsesubarg, 'details')
        tmp <- lapply(tmp, collapsesubarg, 'start')
        tmp <- lapply(tmp, function(x) sapply(x, format))
        nm <- unique(unlist(lapply(tmp, names)))
        tmp0 <- matrix('', nrow = length(tmp), ncol = length(nm),
                       dimnames = list(NULL,  nm))
        ## 2015-05-14 need better labelling of e.g. fixed$sigma, start$g0, start$D
        ## 2016-03-05 ad hoc fix... replace 'core' mask/matrix with scalar
   
        tmp <- lapply(tmp, function(x) {
            if ("core" %in% names(x))
            if (length(x$core)>1)
                x$core <- 'core'
            x})
        for (i in 1:length(tmp)) 
            # tmp0[i,names(tmp[[i]])] <- tmp[[i]]
            tmp0[i,names(tmp[[i]])] <- unlist(tmp[[i]])
        ## as.data.frame(t(unlist(tmp0)))  ## but this fails! 2015-11-03
        as.data.frame(tmp0)
    }
}

header <- function (object) {

    values <- lapply(object$scenarios, unique)
    nvalues <- sapply(values, length)
    notID <- names(object$scenarios) != 'scenario'
    constant <- object$scenarios[1, (nvalues==1) & notID]
    rownames(constant) <- 'value'
    constant <- data.frame(t(constant))
    constant[,1] <- as.character(constant[,1])
    varying <- object$scenarios[, nvalues>1, drop = FALSE]

    ti <- unique(object$scenarios$trapsindex)
    detectors <- data.frame(trapsindex = ti, trapsname = names(object$trapset)[ti])

    pi <- unique(object$scenarios$popindex)
    popargs <- argdf(object$pop.args)
    if (!is.null(popargs))
        popargs <- data.frame(popindex = pi, popargs[pi,,drop = FALSE])

    di <- unique(object$scenarios$detindex)

    detargs <- argdf(object$det.args)
    if (!is.null(detargs))
        detargs <- data.frame(detindex = di, detargs[di,,drop = FALSE])

    fi <- unique(object$scenarios$fitindex)

    fitargs <- argdf(object$fit.args)
    if (!is.null(fitargs))
        fitargs <- data.frame(fitindex = fi, fitargs[fi,,drop = FALSE])

    list(call = object$call, starttime = object$starttime, proctime = object$proctime,
         constant = constant, varying = varying, pop.args = popargs, det.args = detargs,
         fit.args = fitargs, detectors = detectors, nrepl = object$nrepl,
         outputclass = class(object))
}

summary.secrdesign <- function (object, ...) {
    out <- list(header = header(object))
    class (out) <- c('summarysecrdesign', 'list')
    out
 }

summary.rawdata <- function (object, ...) {
    tmp <- fit.models(object)
    summary(tmp, ...)
 }

summary.estimatetables <- function (object, ...) {
    tmp <- select.stats(object)
    summary(tmp, ...)
}

summary.selectedstatistics <- function (object,
                                fields = c('n', 'mean', 'se'),
                                dec = 5,
                                alpha = 0.05,
                                type = c('list','dataframe','array'),
                                ...) {
    if (length(fields) == 1)
        if (tolower(fields) == 'all')
            fields <- c('n', 'mean', 'sd', 'se', 'min', 'max', 'lcl', 'ucl', 'q025',
                        'median', 'q975', 'rms')
    z      <- abs(qnorm(1-alpha/2))   ## beware confusion with hazard z!
    alpha1 <- 0.025
    alpha3 <- 0.975
    qfields <- fields[substring(fields,1,1)=='q']
    if (length(qfields)>0) alpha1 <- as.numeric(substring(qfields[1],2,4))/1000
    if (length(qfields)>1) alpha3 <- as.numeric(substring(qfields[2],2,4))/1000
    q1tag <- paste('q', formatC(1000*alpha1, width=3, format="d", flag = "0"), sep='')
    q3tag <- paste('q', formatC(1000*alpha3, width=3, format="d", flag = "0"), sep='')

    type <- match.arg(type)

    ## allow for zero-length
    minNA <- function (x) if (sum(!is.na(x)) > 0) min(x, na.rm = T) else NA
    maxNA <- function (x) if (sum(!is.na(x)) > 0) max(x, na.rm = T) else NA

    ## data.frame output
    if (type %in% c('list','dataframe')) {
        sumx <- function (x) {
            n    <- sum(!is.na(x))
            mean <- mean(x, na.rm = T)
            sd   <- sd(x, na.rm = T)
            se   <- sd/n^0.5
            minx <- minNA(x)
            maxx <- maxNA(x)
            lcl  <- mean - z * se
            ucl  <- mean + z * se
            q    <- quantile(x, na.rm = T, prob = c(alpha1, 0.5, alpha3))
            rms  <- mean(x^2, na.rm = T)^0.5
            tmp  <- c(n, mean, sd, se, minx, maxx, lcl, ucl, q, rms)
            tmp[is.nan(tmp)] <- NA
            names(tmp) <- c('n', 'mean', 'sd', 'se', 'min', 'max', 'lcl', 'ucl', q1tag,
                            'median', q3tag, 'rms')
            tmp
        }
        valstring <- function (slist) {
            slist <- t(slist)   ## to reorder fields
            val <- as.numeric(slist)
            names(val) <- apply(do.call(expand.grid, dimnames(slist)), 1, paste,
                                collapse='.')
            val
        }
        subscenariolabel <- function (scen) {
            sv <- scen[varying]
            sv <- sapply(sv, format)
            nv <- names(object$scenarios)[varying]
            label <- paste(nv, sv, sep=' = ')
            paste(label, collapse=', ')
        }
        tidy <- function(x) {
            x <- as.data.frame(x)
            x[,fields, drop=F]
        }
        varying <- apply(object$scenarios,2, function(x) length(unique(x)))>1
        varying[1] <- FALSE ## drop scenario code

        tmp <- lapply (object$output, function(y) round(t(apply(y,2,sumx)), dec))
        tmp <- lapply(tmp, tidy)

        if (type == 'list') {
            lab <- apply(object$scenario,1, subscenariolabel)
            names(tmp) <- sapply(split(lab, object$scenarios$scenario), paste,
                                 collapse = ' + ')
        }
        else {
            new <- do.call(rbind, lapply(tmp, valstring))
            tmp <- cbind(object$scenarios[,varying], new)
        }
    
        out <- list(header = header(object), OUTPUT = tmp)
    }
    ## array output
    else {
        outputarray <- make.array(object)
        nd     <- length(dim(outputarray))
        n      <- apply(outputarray, 2:nd, function(y) sum(!is.na(y)))
        mean   <- apply(outputarray, 2:nd, mean, na.rm = T)
        rms    <- apply(outputarray, 2:nd, function (x) mean(x^2, na.rm = T)^0.5)
        sd     <- apply(outputarray, 2:nd, sd, na.rm = T)
        se     <- sd/n^0.5
        minx   <- apply(outputarray, 2:nd, minNA)
        maxx   <- apply(outputarray, 2:nd, maxNA)
        q      <- apply(outputarray, 2:nd, quantile, na.rm = T, prob =
                        c(alpha1, 0.5, alpha3))
        out    <- list(n = n,
                       mean = mean,
                       se = se,
                       sd = sd,
                       min = minx,
                       max = maxx,
                       lcl = mean - z * se,
                       ucl = mean + z * se,
                       rms = rms,
                       q1 = apply(q,2:nd,'[',1),
                       median = apply(q,2:nd,'[',2),
                       q3 = apply(q,2:nd,'[',3))
        names(out)[names(out)=='q1'] <- q1tag
        names(out)[names(out)=='q3'] <- q3tag
        qi <- grepl('q', fields)
        if (any(qi)) {
            fields <- fields[!qi]
            fields <- c(fields, q1tag, q3tag)
        }
        out    <- out[fields]
        statlast <- function (y) {
            y[is.nan(y)] <- NA
            y <- round(y, dec)
            aperm(y, c(2:(nd-1),1))
        }
        tmp    <- lapply(out, statlast)
        out <- list(header = header(object), OUTPUT = tmp)
    }
    ## wrap at 'group'
    out$scenariodetail <- unlist(degroup(names(out$OUTPUT)))
    names(out$OUTPUT) <- as.character(unique(object$scenarios$scenario))
    class (out) <- c('summarysecrdesign', 'list')
    out
}

degroup <- function (x) {
    if (grepl('group', x[[1]])) {
        x <- strsplit(x, 'group')
        x <- lapply(x, function(y) y[nchar(y)>1])
        x <- lapply(x, function(y) paste('group', y, sep=''))
        lapply(x, paste, collapse = '\n ')
    }
    else x
}

print.summarysecrdesign <- function (x, ...) {
    print(x$header$call)
    cat('\n')
    cat('Replicates   ', x$header$nrepl, '\n')
    cat('Started      ', x$header$starttime, '\n')
    cat('Run time     ', round(x$header$proctime/60,3), ' minutes \n')
    cat('Output class ', x$header$outputclass[1], '\n')
    cat('\n')
    print(x$header['constant'])
    cat("$varying\n")
    print(x$header[['varying']], row.names = FALSE)
    cat("\n")
    cat("$detectors\n")
    print(x$header[['detectors']], row.names = FALSE)
    cat("\n")
    if (!is.null(x$header[['pop.args']])) {
        cat("$pop.args\n")
        print (x$header[['pop.args']], row.names = FALSE)
        cat("\n")
    }
    if (!is.null(x$header[['det.args']])) {
        cat("$det.args\n")
        print (x$header[['det.args']], row.names = FALSE)
        cat("\n")
    }
    if (!is.null(x$header[['fit.args']])) {
        cat("$fit.args\n")
        print (x$header[['fit.args']], row.names = FALSE)
        cat("\n")
    }

    if (!is.null(x$OUTPUT)) {
        cat('OUTPUT\n')
        ## replaced 2015-01-27
        ## print(x$OUTPUT)
        for (i in 1:length(x$OUTPUT)) {
            cat ('\n$', names(x$OUTPUT)[i], sep = '')
            if (is.null(x$scenariodetail))
                cat ('\n')
            else
                cat ('\n', x$scenariodetail[i], '\n ')
            print(x$OUTPUT[[i]])
        }
    }
}

#####################################
## plot method for secrdesign object

plot.selectedstatistics <- function (x, scenarios, statistic, type = c('hist','CI'),
                                     refline, xlab = NULL, ...) {
    plothist <- function (mat, stat, scen, ...) {
        if (is.character(stat) & !(stat %in% colnames(mat)))
            stop ("requested statistic not in data")
        else if (is.numeric(stat)) {
            if (stat > dim(mat)[2])
                stop ("statistic exceeds dimension of data")
            stat <- colnames(mat)[stat]
        }
        if (is.null(xlab)) {
            if (!is.null(param))
                xlab <- paste(stat, ' (', param, ')', sep = '')
            else
                xlab <- stat
        }
        hist (mat[,stat], main = "", xlab = xlab, ...)
        if (refline) {
            if (param %in% c('E.N','R.N')) {
                true <- x$scenarios[scen, 'D']
                ## true <- true * attr(x, 'regionarea')[scen]
                ## 2014-11-12
                true <- true * attr(x, 'regionsize')[scen]
            }
            else {
                true <- x$scenarios[scen, param]
            }
            abline(v = true, col = 'red')
        }
        mtext(side=3, line = 1, paste('Scenario', scen, collapse=' '), cex = par()$cex*1.1)
    }
    plotCI <- function (mat, scen, ...) {
        estname <- 'estimate'
        if (x$outputtype == 'coef')
            estname <- 'beta'
        if (!all(c(estname,'lcl','ucl') %in% colnames(mat)))
            stop ("requires statistics 'estimate','lcl','ucl'")
        if (is.null(xlab))
            xlab <- 'Replicate'
        if (!is.null(param))
            ylab <- paste(estname, ' (', param, ')', sep='')
        nrepl <- nrow(mat)
        plot (c(1,nrepl), range(mat[,c('lcl','ucl')], na.rm = TRUE), type = 'n',
              xlab = xlab, ylab = ylab, ...)
        segments (1:nrepl, mat[,'lcl'], 1:nrepl, mat[,'ucl'])
        if (refline) {
            if (param %in% c('E.N','R.N')) {
                true <- x$scenarios[scen, 'D']
                ## true <- true * attr(x, 'regionarea')[scen]
                ## 2014-11-12
                true <- true * attr(x, 'regionsize')[scen]
            }
            else {
                true <- x$scenarios[scen, param]
            }
            abline(h = true, col = 'red')
        }
        points (1:nrepl, mat[,estname], pch = 21, bg = 'white')
        mtext(side=3, line = 1, paste('Scenario', scen, collapse=' '), cex = par()$cex*1.1)
    }
    param <- attr(x, 'parameter')
    type <- match.arg(type)
    if (missing(scenarios))
        scenarios <- 1:length(x$output)
    if (missing(statistic))
        statistic <- 1
    if (missing(refline))
        refline <- type == 'CI'   ## default TRUE if CI, FALSE otherwise
    for (scenario in scenarios) {
        if (type == 'CI')
            plotCI(x$output[[scenario]], scen = scenario,...)
        else {
            for (s in statistic) {
                if (type == 'hist')
                    plothist(x$output[[scenario]], stat = s, scen = scenario,...)
            }
        }
    }
}
