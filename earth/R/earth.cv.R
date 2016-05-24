# earth.cv.R: Functions for cross validation of earth models.
#             Note that earth.cv returns null unless nfold > 1.

earth.cv <- function(object, x, y,
    subset, weights, na.action, pmethod, keepxy, trace, glm, degree, nprune,
    ncross, nfold, stratify, get.oof.fit.tab, get.oof.rsq.per.subset,
    Scale.y, env, ...)
{
    get.fold.rsq.per.subset <- function(foldmod, oof.y, max.nterms, trace, must.print.dots)
    {
        wp.expanded <- wp.expanded / sum(wp.expanded)
        oof.rsq.per.subset <- infold.rsq.per.subset <- repl(0, max.nterms)
        # nrow(foldmod$dirs) is the number of terms in this fold's model before pruning
        for(nterms in 1:min(max.nterms, nrow(foldmod$dirs))) {
            trace.get.fold1(trace, must.print.dots, nterms)
            # penalty=-1 to enforce strict nprune TODO consider changing to pmethod=none
            # glm=NULL for speed, ok because we don't need the glm submodel
            # TODO with keepxy=TRUE, 70% of cv time is spent in update.earth
            pruned.foldmod <- update(foldmod, nprune=min(nprune, nterms),
                    penalty=-1, ponly=TRUE, glm=NULL, trace=max(0, trace-1))
            fit <- predict(pruned.foldmod, newdata=x, type="earth")
            oof.fit    <- fit[oof.subset,  , drop=FALSE]
            infold.fit <- fit[infold.subset, , drop=FALSE]
            for(iresp in seq_len(NCOL(fit))) { # for each response
                oof.rsq.per.subset[nterms] <- oof.rsq.per.subset[nterms] +
                    wp.expanded[iresp] *
                        get.weighted.rsq(oof.y[,iresp], oof.fit[,iresp], oof.weights)

                infold.rsq.per.subset[nterms] <- infold.rsq.per.subset[nterms] +
                    wp.expanded[iresp] *
                        get.weighted.rsq(infold.y[,iresp], infold.fit[,iresp], infold.weights)
            }
            if(nrow(foldmod$dirs) < max.nterms)
                for(nterms in (nrow(foldmod$dirs)+1): max.nterms)
                    oof.rsq.per.subset[nterms] <- infold.rsq.per.subset[nterms] <- NA
        }
        trace.get.fold2(trace, must.print.dots, nterms)
        list(oof.rsq.per.subset    = oof.rsq.per.subset,
             infold.rsq.per.subset = infold.rsq.per.subset)
    }
    #--- earth.cv starts here ---
    # We called check.cv.args(ncross, nfold, varmod.method, pmethod)
    # earlier so it's safe to use those args here
    # Likewise, subset arg was already checked in earth.fit.
    stratify <- check.boolean(stratify)
    stopifnot(ncross >= 1, nfold > 1)
    trace1(trace, "\n")
    if(nfold > nrow(x))
        nfold <- nrow(x)
    if(!is.null(object$glm.bpairs))
        stop0("earth does not yet support cross validation of paired binomial responses")
    max.nterms <- nrow(object$dirs)
    wp <- wp.expanded <- object$wp
    if(is.null(wp))
        wp.expanded <- repl(1, ncol(y)) # all ones vector
    cv.list <- list()                   # returned list of cross validated models
    ncases <- nrow(x)
    nresp <- ncol(y)                    # number of responses
    # ndigits aligns trace prints without too much white space
    ndigits <- ceiling(log10(ncases - ncases/nfold + .1))
    # print pacifier dots if get.fold.rss.per.subset will be slow
    must.print.dots <- trace >= .5 && trace <= 1 &&
        (get.oof.rsq.per.subset || get.oof.fit.tab) &&
        nrow(x) * max.nterms > 50e3
    trace.cv.header(object, nresp, trace, must.print.dots)
    n.oof.digits <- ceiling(log10(1.2 * ncases / nfold)) # 1.2 allows for diff sized subsets
    groups <- matrix(NA, nrow=ncross*ncases, ncol=2)
    colnames(groups) <- c("cross", "fold")
    for(icross in seq_len(ncross)) {
        start <- ((icross-1) * ncases) + 1
        groups[start:(start+ncases-1), 1] <- icross
        groups[start:(start+ncases-1), 2] <- get.groups(y, nfold, stratify)
    }
    is.binomial <- is.poisson <- FALSE
    if(!is.null(glm)) {
        glm1 <- get.glm.arg(glm)
        family <- get.glm.family(glm1$family, env=env)
        is.binomial <- is.binomial(family)
        is.poisson <- is.poisson(family)
    }
    must.get.class.rate <- !is.null(object$levels)
    # the final summary row of the tables is "mean", "all" or "max", depending on the statistic.
    if(ncross > 1)
        fold.names <- paste0(rep(paste0("fold", seq_len(ncross), "."), each=nfold),
                             rep(seq_len(nfold), times=ncross))
    else
        fold.names <- sprintf("fold%d", seq_len(nfold))
    fold.names.plus.mean <- c(fold.names, "mean")
    fold.names.plus.all  <- c(fold.names, "all")
    fold.names.plus.max  <- c(fold.names, "max")
    resp.names.plus.mean <- c(colnames(y), "mean") # response names plus "mean"
    resp.names.plus.max  <- c(colnames(y), "max")
    ncross.fold <- ncross * nfold
    nvars.selected.by.gcv <- double(ncross.fold+1)  # nbr of used predictors in each CV mod
    nterms.selected.by.gcv <- double(ncross.fold+1) # nbr of selected terms in each CV mod
    names(nvars.selected.by.gcv) <- fold.names.plus.mean
    names(nterms.selected.by.gcv) <- fold.names.plus.mean

    rsq.tab <- matrix(0, nrow=ncross.fold+1, ncol=1+nresp) # table of cv results, +1 for means
    colnames(rsq.tab) <- resp.names.plus.mean
    rownames(rsq.tab) <- fold.names.plus.mean

    maxerr.tab <- matrix(0, nrow=ncross.fold+1, ncol=1+nresp) # table of cv results, +1 for max
    colnames(maxerr.tab) <- resp.names.plus.max
    rownames(maxerr.tab) <- fold.names.plus.max

    deviance.tab <- calib.int.tab <- calib.slope.tab <- test.tab <- class.rate.tab <- NULL
    if(is.binomial || is.poisson) {
        deviance.tab    <- matrix(0, nrow=ncross.fold+1, ncol=1+nresp) # mean deviance
        calib.int.tab   <- matrix(0, nrow=ncross.fold+1, ncol=1+nresp)
        calib.slope.tab <- matrix(0, nrow=ncross.fold+1, ncol=1+nresp)
        test.tab        <- matrix(0, nrow=ncross.fold+1, ncol=1+nresp) # binomial auc, poisson cor
        colnames(deviance.tab) <- colnames(calib.int.tab) <-
            colnames(calib.slope.tab) <- colnames(test.tab) <- resp.names.plus.mean
        rownames(deviance.tab) <- rownames(calib.int.tab) <-
            rownames(calib.slope.tab) <- rownames(test.tab) <- fold.names.plus.mean
    }
    if(must.get.class.rate) {
        class.rate.tab <- matrix(0, nrow=ncross.fold+1, ncol=1+nresp)
        colnames(class.rate.tab) <- resp.names.plus.mean
        rownames(maxerr.tab) <- fold.names.plus.all
    }
    oof.rsq.tab <- infold.rsq.tab <- oof.fit.tab <- NULL
    if(get.oof.rsq.per.subset) {
        oof.rsq.tab <- infold.rsq.tab <- matrix(0, nrow=ncross.fold+1, ncol=max.nterms)
        colnames(oof.rsq.tab) <- colnames(infold.rsq.tab) <- paste0("nterms", 1:max.nterms)
        rownames(oof.rsq.tab) <- rownames(infold.rsq.tab) <- fold.names.plus.all
    }
    if(get.oof.fit.tab) {
        oof.fit.tab <- matrix(0, nrow=nrow(x), ncol=ncross) # preds on oof data
        colnames(oof.fit.tab) <- paste0("icross", seq_len(ncross))
    }
    for(icross in seq_len(ncross)) {
        this.group <- get.this.group(icross, ifold, ncases, groups)
        for(ifold in seq_len(nfold)) {
            icross.fold <- ((icross-1) * nfold) + ifold
            oof.subset <- seq_len(ncases)[which(this.group == ifold)]
            infold.subset <- seq_len(ncases)[-oof.subset]
            trace.fold.header(trace, ncross, icross, ifold)
            infold.x <- x[infold.subset,,drop=FALSE]
            infold.y <- y[infold.subset,,drop=FALSE]
            infold.weights <- weights[infold.subset]
            foldmod <- earth.default(x=infold.x, y=infold.y, weights=infold.weights,
                wp=wp, Scale.y=Scale.y, subset=subset,
                trace=trace, glm=glm, degree=degree,
                pmethod=if(pmethod == "cv") "backward" else pmethod,
                ncross=0, nfold=0, varmod.method="none",
                ...)
            foldmod$icross <- icross
            foldmod$ifold <- ifold
            oof.x <- x[oof.subset,,drop=FALSE]
            oof.y <- y[oof.subset,,drop=FALSE]
            oof.weights <- weights[oof.subset]
            oof.fit <- predict(foldmod, newdata=oof.x, type="earth")
            # fill in subset of entries in this icross column of oof.fit.tab
            # note that we use only the first response when there are multiple responses
            if(!is.null(oof.fit.tab))
                oof.fit.tab[oof.subset, icross] <- oof.fit[,1]
            oof.fit.resp <- NULL
            if(is.binomial || is.poisson)
                oof.fit.resp <- predict(foldmod, newdata=oof.x, type="response")
            else if(must.get.class.rate) # not glm but has binary response?
                oof.fit.resp <- oof.fit

            # fill in this fold's row in summary tabs

            for(iresp in seq_len(nresp)) {
                rsq.tab[icross.fold, iresp] <-
                    get.weighted.rsq(oof.y[,iresp], oof.fit[,iresp], oof.weights)
                if(is.binomial) {
                    deviance.tab[icross.fold, iresp] <-
                        get.binomial.deviance(oof.fit.resp[,iresp], oof.y[,iresp])
                    calib <- get.binomial.calib(oof.fit.resp[,iresp], oof.y[,iresp])
                    calib.int.tab[icross.fold, iresp]   <- calib[1]
                    calib.slope.tab[icross.fold, iresp] <- calib[2]
                    maxerr.tab[icross.fold, iresp] <-
                        get.maxerr(oof.y[,iresp] - oof.fit.resp[,iresp])
                    test.tab[icross.fold, iresp] <-
                        get.auc(oof.fit.resp[,iresp], oof.y[,iresp])
                }
                else if(is.poisson) {
                    deviance.tab[icross.fold, iresp] <-
                        get.poisson.deviance(oof.fit.resp[,iresp], oof.y[,iresp])
                    calib <- get.poisson.calib(oof.fit.resp[,iresp], oof.y[,iresp])
                    calib.int.tab[icross.fold, iresp]   <- calib[1]
                    calib.slope.tab[icross.fold, iresp] <- calib[2]
                    maxerr.tab[icross.fold, iresp] <-
                        get.maxerr(oof.y[,iresp] - oof.fit.resp[,iresp])
                    test.tab[icross.fold, iresp] <-
                        cor(oof.fit.resp[,iresp], oof.y[,iresp])
                }
                else
                    maxerr.tab[icross.fold, iresp] <-
                        get.maxerr(oof.y[,iresp] - oof.fit[,iresp])
            } # end for iresp

            nvars.selected.by.gcv[icross.fold] <-
                get.nused.preds.per.subset(foldmod$dirs, foldmod$selected.terms)
            nterms.selected.by.gcv[icross.fold] <-
                length(foldmod$selected.terms)
            if(must.get.class.rate)
                class.rate.tab[icross.fold, ] <- get.class.rate(oof.fit.resp, oof.y, object$levels)
            if(get.oof.rsq.per.subset) {
                temp <-
                    get.fold.rsq.per.subset(foldmod, oof.y, max.nterms, trace, must.print.dots)
                oof.rsq.tab[icross.fold,]    <- temp$oof.rsq
                infold.rsq.tab[icross.fold,] <- temp$infold.rsq
            }
            # init last column of summary tables
            ilast.col <- nresp+1 # index of final (summary) column in tables
            rsq.tab[icross.fold, ilast.col]    <- weighted.mean(rsq.tab[icross.fold, -ilast.col], wp.expanded)
            maxerr.tab[icross.fold, ilast.col] <- get.maxerr(maxerr.tab[icross.fold, -ilast.col])
            if(is.binomial || is.poisson) {
                deviance.tab[icross.fold, ilast.col]    <- weighted.mean(deviance.tab   [icross.fold, -ilast.col], wp.expanded)
                calib.int.tab[icross.fold, ilast.col]   <- weighted.mean(calib.int.tab  [icross.fold, -ilast.col], wp.expanded)
                calib.slope.tab[icross.fold, ilast.col] <- weighted.mean(calib.slope.tab[icross.fold, -ilast.col], wp.expanded)
                test.tab[icross.fold, ilast.col]        <- weighted.mean(test.tab       [icross.fold, -ilast.col], wp.expanded)
            }
            if(!keepxy) # reduce memory by getting rid of big fields
                foldmod$bx <- foldmod$residuals <- foldmod$prune.terms <- NULL
            trace.fold(icross, ifold, trace, y, infold.subset, oof.subset, ncross, ndigits,
                       rsq.tab[icross.fold,], n.oof.digits, must.print.dots)
            cv.list[[icross.fold]] <- foldmod
        } # end for ifold
    } # end for icross

    # init last row of summary tables

    ilast <- ncross.fold+1 # index of last row in tables
    nvars.selected.by.gcv [ilast]  <- mean(nvars.selected.by.gcv[-ilast])
    nterms.selected.by.gcv[ilast]  <- mean(nterms.selected.by.gcv[-ilast])

    rsq.tab   [ilast,] <- colMeans(rsq.tab     [-ilast,])
    maxerr.tab[ilast,] <- get.maxerr(maxerr.tab[-ilast,])
    if(is.binomial || is.poisson) {
        deviance.tab   [ilast,] <- colMeans(deviance.tab   [-ilast,])
        calib.int.tab  [ilast,] <- colMeans(calib.int.tab  [-ilast,])
        calib.slope.tab[ilast,] <- colMeans(calib.slope.tab[-ilast,])
        test.tab       [ilast,] <- colMeans(test.tab       [-ilast,])
    }
    if(must.get.class.rate)
        class.rate.tab[ilast,]  <- colMeans(class.rate.tab [-ilast,])

    oof.rsq.per.subset <- NULL
    if(get.oof.rsq.per.subset) {
        # there will be NAs in oof.rsq.tab if max terms in a fold is
        # less than max terms in full model
        oof.rsq.tab[ilast,]    <-
            col.means.with.special.na.handling(oof.rsq.tab[-ilast,])
        infold.rsq.tab[ilast,] <-
            col.means.with.special.na.handling(infold.rsq.tab[-ilast,])
    }
    trace1(trace, "\n")
    trace.fold(icross, -1, trace, y, TRUE, TRUE, ncross, ndigits,
               rsq.tab[ilast,], n.oof.digits, must.print.dots)
    if(trace >= .5)
        cat("\n")
    names(cv.list) <- fold.names
    rv <- list(
        cv.list = cv.list, # list of earth models built during cross validation
        nterms.selected.by.gcv = nterms.selected.by.gcv,
        nvars.selected.by.gcv  = nvars.selected.by.gcv,
        groups          = groups, # groups used for cross validation
        rsq.tab         = rsq.tab,
        maxerr.tab      = maxerr.tab,
        class.rate.tab  = class.rate.tab,
        auc.tab         = if(is.binomial) test.tab else NULL,
        cor.tab         = if(is.poisson)  test.tab else NULL,
        deviance.tab    = deviance.tab,
        calib.int.tab   = calib.int.tab,
        calib.slope.tab = calib.slope.tab,
        oof.fit.tab     = oof.fit.tab,
        infold.rsq.tab  = infold.rsq.tab,
        oof.rsq.tab     = oof.rsq.tab)
    rv
}
# Return the mean of each column in tab.
# NAs are ignored, except that the column mean is NA
# for columns in which over half the entries are NA.

col.means.with.special.na.handling <- function(tab)
{
    means <- colMeans(tab, na.rm=TRUE)
    nna <- colSums(is.na(tab)) # number of NAs in each column
    means[nna > nrow(tab) / 2] <- NA
    means # a vector of column means, some entries may be NA
}
check.cv.args <- function(ncross, nfold, pmethod, varmod.method)
{
    check.integer.scalar(ncross, min=0)
    if(ncross > 1000) # 1000 is arbitrary
        stop0("ncross ", ncross, " is too big")

    check.integer.scalar(nfold, min=0)
    # if(nfold > 10000) # 10000 is arbitrary
    #     stop0("nfold ", nfold, " is too big")

    if(ncross > 1 && nfold < 2)
        stop0("ncross=", ncross, " yet nfold=", nfold)
    if(ncross < 1 && nfold > 1)
        stop0("ncross=", ncross, " yet nfold=", nfold)

    if(varmod.method != "none") {
        if(nfold <= 1)
            stop0("varmod.method=\"", varmod.method, "\" requires nfold greater than 1")
        if(ncross < 3)
            stop0("ncross=", ncross,
                  " but should be larger when varmod.method is used\n",
                  "       (suggest at least 30, for debugging 3 is ok)")
    }
}
get.groups <- function(y, nfold, stratify)
{
    groups <- sample(repl(seq_len(nfold), nrow(y)))
    if(stratify) {
        # Get (roughly) equal number of folds for each non-zero entry in each y column
        # If y was originally a factor before expansion to multiple columns, this is
        # equivalent to having the same numbers of each factor level in each fold.
        for(iresp in seq_len(ncol(y))) {
            yset <- y[,iresp] != 0
            groups[yset] <- sample(repl(seq_len(nfold), sum(yset)))
        }
    }
    if(any(table(groups) == 0))
        stop0("Not enough data to do ", nfold,
              " fold cross validation (an out-of-fold set is empty)")
    groups
}
get.this.group <- function(icross, ifold, ncases, groups)
{
    start <- ((icross-1) * ncases) + 1
    groups[start:(start+ncases-1), 2]
}
trace.cv.header <- function(object, nresp, trace, must.print.dots)
{
    if(trace == .5) {
        if(must.print.dots || nresp > 1)
            cat("\n")
        printf("Full model GRSq %5.3f RSq %5.3f, starting cross validation\n",
               object$grsq , object$rsq)
        if(must.print.dots || nresp > 1)
            cat("\n")
    }
}
trace.fold.header <- function(trace, ncross, icross, ifold)
{
    if(trace >= .5 && trace < 1) {
        if(ncross > 1)
            printf("CV fold %2d.%-2d ", icross, ifold)
        else
            printf("CV fold %-2d ", ifold)
    } else if(trace >= 1) { # newline etc. to distinguish this from other trace prints
        if(ncross > 1)
            printf("\nCV fold %d.%d -----------------------------%s", icross, ifold,
                   "---------------------------------------\n")
        else
            printf("\nCV fold %d -----------------------------%s", ifold,
                   "---------------------------------------\n")
    }
}
# print results for the current fold (ifold=-1 means "all")

trace.fold <- function(icross, ifold, trace, y, infold.subset, oof.subset, ncross, ndigits,
                       rsq.row, n.oof.digits, must.print.dots)
{
    if(trace < .5)
        return()
    if(ifold < 0) {
        icross <- if(ncross > 1) "   " else ""
        printf("%s%s", "CV all     ", icross)
    } else if(trace >= 1) {
        if(ncross > 1)
            icross <- sprintf("%2d.", icross)
        else
            icross <- ""
        printf("CV fold %s%-2d ", icross, ifold)
    }
    printf("CVRSq %-6.3f ", rsq.row[length(rsq.row)])
    nresp <- length(rsq.row) - 1 # -1 for "all"
    if(nresp > 1) {
        cat("Per response CVRSq ")
        for(iresp in seq_len(nresp))
            printf("%-6.3f ", rsq.row[iresp])
    }
    if(nresp > 1) {
        if(ncross > 1)
            cat("\n                            ")
        else
            cat("\n                         ")
    }
    if(ifold < 0)
        printf("      %.*s       ", n.oof.digits, "                    ")
    else
        printf("n.oof %*.0f %2.0f%%  ", n.oof.digits,
               length(infold.subset),
               100 * (length(y[, 1]) - length(infold.subset)) / length(y[, 1]))
    cat("n.infold.nz ")
    if(nresp == 1)
        printf("%*.0f %2.0f%%  ", ndigits, sum(y[infold.subset, 1] != 0),
               100 * sum(y[infold.subset, 1] != 0) / length(y[infold.subset, 1]))
    else for(iresp in seq_len(nresp))
        printf("%*.0f ", ndigits, sum(y[infold.subset, iresp] != 0))

    if(ifold >= 0) {
        cat("n.oof.nz ")
        if(nresp == 1)
            printf("%*.0f %2.0f%%  ", ndigits, sum(y[oof.subset, 1] != 0),
                   100 * sum(y[oof.subset, 1] != 0) / length(y[oof.subset, 1]))
        else for(iresp in seq_len(nresp))
            printf("%*.0f  ", ndigits, sum(y[oof.subset, iresp] != 0))
        if(nresp > 1)
            cat("\n")
    }
    if(must.print.dots && trace == .5 && ifold > 0)
        cat("\n")
    cat("\n")
}
trace.get.fold1 <- function(trace, must.print.dots, nterms)
{
    if(trace >= 2)
        cat0("\nget.fold.rss.per.subset nterms=", nterms, "\n")
    else if(must.print.dots) {
        cat0(".")
        if(nterms %% 40 == 0) {
            cat0("\n")
            if(trace == .5)
                cat("            ")
        }
        flush.console()
    }
}
trace.get.fold2 <- function(trace, must.print.dots, nterms)
{
    if(must.print.dots && nterms %% 40) { # nterms %% 40 avoids double newline
        cat0("\n")
        if(trace == .5)
            cat("            ")
    }
}
print.cv <- function(x) # called from print.earth for cross validated models
{
    cv.field <- function(field) round(x[[field]][ilast,], 3)

    cv.sd <- function(field) round(apply(x[[field]][-ilast,], 2, sd), 3)

    #--- print.cv starts here ---
    stopifnot(!is.null(x$cv.list), x$pmethod != "cv")
    ilast <- nrow(x$cv.rsq.tab) # index of "all" row in summary tables
    cat("\nNote: the cross-validation sd's below are standard deviations across folds\n\n")

    printf("Cross validation:   nterms %.2f sd %.2f    nvars %.2f sd %.2f\n\n",
        x$cv.nterms.selected.by.gcv[ilast],
        sd(x$cv.nterms.selected.by.gcv[-ilast]),
        x$cv.nvars.selected.by.gcv[ilast],
        sd(x$cv.nvars.selected.by.gcv[-ilast]))

    # if printing little then use wide spacing
    wide.spacing <- is.null(x$cv.deviance.tab)

    # create a data.frame and print that
    tab <- if(wide.spacing)
               data.frame("    CVRSq"=cv.field("cv.rsq.tab"), check.names = FALSE)
           else
               data.frame("CVRSq"  =cv.field("cv.rsq.tab"))

    tab$sd.1=cv.sd("cv.rsq.tab")

    if(!is.null(x$cv.class.rate.tab)) {
        if(wide.spacing) tab$"    ClassRate" <- cv.field("cv.class.rate.tab")
        else             tab$"ClassRate"   <- cv.field("cv.class.rate.tab")

        tab$sd.2 <- cv.sd("cv.class.rate.tab")
    }
    if(wide.spacing)
        tab$"    MaxErr" <- x$cv.maxerr.tab[ilast,]
    else
        tab$"MaxErr"   <- x$cv.maxerr.tab[ilast,]

    tab$sd <- apply(x$cv.maxerr.tab[-ilast,], 2, sd)

    if(!is.null(x$cv.auc.tab)) {
        tab$AUC  <- cv.field("cv.auc.tab")
        tab$sd.3 <- cv.sd("cv.auc.tab")
    }
    if(!is.null(x$cv.cor.tab)) {
        tab$cor.tab <- cv.field("cv.cor.tab")
        tab$sd.4    <- cv.sd("cv.cor.tab")
    }
    if(!is.null(x$cv.deviance.tab)) {
        tab$MeanDev    <- x$cv.deviance.tab[ilast,]
        tab$sd.5       <- apply(x$cv.deviance.tab[-ilast,], 2, sd)
        tab$CalibInt   <- cv.field("cv.calib.int.tab")
        tab$sd.6       <- cv.sd("cv.calib.int.tab")
        tab$CalibSlope <- cv.field("cv.calib.slope.tab")
        tab$sd.7       <- cv.sd("cv.calib.slope.tab")
    }
    # change "sd.N" to plain "sd"
    names <- names(tab)
    names[grep("sd", names)] <- "sd"
    names(tab) <- names

    rownames(tab) <- c(colnames(x$fitted.values), "All")

    digits <- min(getOption("digits"), 3)

    if(NCOL(x$coefficients) == 1)   # single response model?
        print(tab[1,,drop=FALSE], digits=digits, row.names=FALSE) # skip "All" row
    else                            # multiple response model
        print(tab, digits=digits)
}
