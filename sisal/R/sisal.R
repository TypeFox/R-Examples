### File R/sisal.R
### This file is part of the sisal package for R.
###
### Copyright (C) 2015 Aalto University
###
### This program is free software; you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation; either version 2 of the License, or
### (at your option) any later version.
###
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
### GNU General Public License for more details.
###
### A copy of the GNU General Public License is available at
### http://www.r-project.org/Licenses/

sisal <- function(X, y, Mtimes=100, kfold=10, hbranches=1,
                  max.width=hbranches^2, q=0.165, standardize=TRUE,
                  pruning.criterion=c("round robin", "random nodes",
                      "random edges", "greedy"),
                  pruning.keep.best=TRUE, pruning.reverse=FALSE,
                  verbose=1, use.ridge=FALSE,
                  max.warn=getOption("nwarnings"), sp=-1, ...) {

    pruning.criterion2 <- match.arg(pruning.criterion)
    prune.rr <- pruning.criterion2 == "round robin"
    prune.edges <- pruning.criterion2 == "random edges"
    prune.greedy <- pruning.criterion2 == "greedy"
    stopifnot(identical(pruning.keep.best, TRUE) ||
              identical(pruning.keep.best, FALSE),
              identical(pruning.reverse, TRUE) ||
              identical(pruning.reverse, FALSE),
              identical(standardize, TRUE) || identical(standardize, FALSE),
              identical(use.ridge, TRUE) || identical(use.ridge, FALSE),
              is.numeric(verbose), length(verbose) == 1,
              is.numeric(hbranches), length(hbranches) == 1,
              round(hbranches) == hbranches, hbranches >= 1,
              is.numeric(max.width), length(max.width) == 1,
              round(max.width) == max.width, max.width >= 1,
              is.numeric(q), length(q) == 1, q > 0, q < 0.5,
              is.numeric(Mtimes), length(Mtimes) == 1,
              round(Mtimes) == Mtimes, Mtimes >= 1,
              is.numeric(max.warn), length(max.warn) == 1,
              round(max.warn) == max.warn, max.warn >= 0,
              is.numeric(sp), length(sp) == 1, is.finite(sp))
    allWarn <- vector(mode="list", length=max.warn)
    n.warn <- 0
    X2 <- as.matrix(X)
    if (!is.matrix(X2)) {
        stop(gettextf("'%s' must be coercible to a matrix", "X",
                      domain = "R-sisal"), domain = NA)
    }
    if (!is.numeric(X2)) {
        stop(gettextf("'%s' must be numeric", "X", domain = "R-sisal"),
             domain = NA)
    }
    n <- nrow(X2)
    d <- ncol(X2)
    if (d == 0) {
        stop(gettextf("'%s' must have at least one column", "X",
                      domain = "R-sisal"), domain = NA)
    }
    y2 <- as.vector(y)
    if (!is.vector(y2) || length(y2) != n) {
        stop(gettextf("'%s' must be coercible to a vector with 'length(%s) == nrow(%s)'",
                      "y", "y", "X", domain = "R-sisal"), domain = NA)
    }
    if (!is.numeric(y2)) {
        stop(gettextf("'%s' must be numeric", "y", domain = "R-sisal"),
             domain = NA)
    }

    orig.names <- colnames(X)
    if (is.null(orig.names)) {
        X.names <- paste0("V", 1:d)
    } else {
        X.names <- orig.names
    }
    mean.y <- mean(y2, na.rm=TRUE)
    sd.y <- sd(y2, na.rm=TRUE)
    mean.X <- .colMeans(X2, n, d, na.rm=TRUE)
    sd.X <- apply(X2, 2, sd, na.rm=TRUE)
    if ("(Intercept)" %in% orig.names) {
        stop('"(Intercept)" is a forbidden variable name')
    }
    names(mean.X) <- X.names
    names(sd.X) <- X.names
    zeroRange.X <- logical(d)
    for (k in seq_len(d)) {
        zeroRange.X[k] <- zeroRange(X2[, k], mean.x = mean.X[k])
    }
    zeroRange.y <- zeroRange(y2, mean.x = mean.y)

    isZero <- function (x) {
        isTRUE(all.equal(0, x, check.attributes = FALSE))
    }

    ## Columns of X that are marked zero or constant may still have NA values
    zero.X <- zeroRange.X
    zero.X[zero.X] <- vapply(mean.X[zero.X], isZero, TRUE)
    xNA <- is.na(X2)
    constant.X <- zeroRange.X & !zero.X & colSums(!xNA) > 0
    n.zeros <- sum(zero.X)
    n.constants <- sum(constant.X)

    add.constant <- n.constants == 0
    if (!add.constant) {
        message(sprintf(ngettext(n.constants,
                                 "%.0f constant column found in 'X', will take part in input selection\n",
                                 "%.0f constant columns found in 'X', will take part in input selection\n",
                                 domain = "R-sisal"),
                        n.constants),
                domain = NA)
    }
    if (n.constants > 1) {
        tmp <-
            quote(gettext("selection between >1 constant columns is arbitrary",
                          domain = "R-sisal"))
        n.warn <- n.warn + 1
        if (n.warn <= max.warn) {
            allWarn[[n.warn]] <- tmp
        }
        warning(eval(tmp), domain = NA)
    }

    if (standardize) {
        y2 <- y2 - mean.y
        if (!zeroRange.y && sd.y != 0) {
            y2 <- y2 / sd.y
        }
        sd.X.tmp <- sd.X
        sd.X.tmp[zeroRange.X | sd.X == 0] <- 1
        mean.X.tmp <- mean.X
        mean.X.tmp[constant.X] <- mean.X[constant.X] - 1
        X2 <- sweep(sweep(X2, 2, mean.X.tmp, "-"), 2, sd.X.tmp, "/")
    } else if (use.ridge &&
               (!isTRUE(all.equal(rep.int(sd.X[1], d), sd.X,
                                  check.attributes = FALSE)) ||
                !isTRUE(all.equal(rep.int(0, d + 1), c(mean.X, mean.y),
                                  check.attributes = FALSE)))) {
        warning("the data should be standardized when 'use.ridge' is TRUE")
    }

    ## X.complete and y.complete are possibly standardized versions of X
    ## (with column names) and y.  All samples are still in place.
    X.complete <- X2
    colnames(X.complete) <- X.names
    y.complete <- y2

    ## Remove NA samples (y or all columns of X are NA)
    sumNA <- rowSums(xNA)
    naX <- sumNA == d
    naY <- is.na(y2)
    naRow <- naX | naY
    n.missing <- sum(naRow)
    if (n.missing > 0) {
        if (verbose > 0) {
            splitcat(gettextf("%d samples dropped due to missing values",
                              n.missing, domain = "R-sisal"))
        }
        n <- n - n.missing
        notNA <- !naRow
        X2 <- X2[notNA, , drop = FALSE]
        y2 <- y2[notNA]
        xNA <- xNA[notNA, , drop = FALSE]
        sumNA <- sumNA[notNA]
    }

    goodX0 <- rowSums(xNA) == 0
    notGood0 <- which(!goodX0)
    nGood0 <- sum(goodX0)
    idx.clean <- sumNA == 0
    n.clean <- sum(idx.clean)
    ## 2 complete samples is, of course, not really enough. However,
    ## an assumption made later in the code is now fulfilled: at least
    ## one set of model coefficients (even if poorly estimated) will
    ## always be available.
    if (n.clean < 2) {
        stop("too few complete samples: 2 is the absolute minimum")
    }

    ## Elementary sanity checks for values of Mtimes and kfold
    if (!is.numeric(kfold) || length(kfold) != 1 || round(kfold) != kfold ||
        kfold > n || kfold < 2) {
        stop("'kfold' must be an integral valued number, 2 <= 'kfold' <= number of samples")
    }
    if(kfold == n && Mtimes > 1) {
        stop("Leave-one-out cross-validation should not be repeated.",
             "Adjust 'Mtimes' and/or 'kfold'.")
    } else if (kfold == n - 1 && Mtimes > n * (n - 1) / 2) {
        stop("At least two cross-validation repetitions would be identical. ",
             "Adjust 'Mtimes' and/or 'kfold'.")
    } else if (kfold == n - 2 &&
               Mtimes > n * (n - 1) * (n - 2) * (n - 3) / 8) {
        stop("At least two cross-validation repetitions would be identical. ",
             "Adjust 'Mtimes' and/or 'kfold'.")
    }

    params <- c(mget(c("pruning.criterion2",
                       setdiff(names(formals()),
                               c("X", "y", "verbose", "max.warn",
                                 "pruning.criterion", "..."))),
                     environment()),
                list(...))
    names(params)[names(params) == "pruning.criterion2"] <-
        "pruning.criterion"
    if (verbose > 0) {
        splitcat(paste(names(params), params, sep="=", collapse=", "))
        splitcat(gettextf("There are %.0f variables and %.0f samples.",
                          d, n, domain = "R-sisal"))
        if (n.clean == n && n.missing == 0) {
            splitcat(gettext("There are no missing values.",
                             domain = "R-sisal"))
        } else if (n.clean == n) {
            splitcat(gettext("There are no missing values in the remaining samples.",
                             domain = "R-sisal"))
        } else {
            n.dirty <- n - n.clean
            splitcat(sprintf(ngettext(n.dirty,
                                      "There is %.0f sample with one or more missing values.",
                                      "There are %.0f samples with one or more missing values.",
                                      domain = "R-sisal"),
                             n.dirty))
        }
    }

    floor.size <- floor(n / kfold)
    part.sizes <- rep.int(floor.size, kfold)
    part.sizes[seq_len(n - kfold * floor.size)] <- floor.size + 1
    part.end <- cumsum(part.sizes)
    part.start <- part.end - part.sizes + 1

    selected <- rep.int(TRUE, d)
    n.selected <- d
    D <- d + 1                          # + 1 to account for the empty model
    E.v.best <- rep.int(0, D)
    E.tr.best <- rep.int(0, D)
    s.tr.best <- rep.int(0, D)
    names.best <- character(D)

    rank.required <- d + 1 - n.constants - n.zeros
    min.subproblem <- n - max(part.sizes)
    ## Quick test.  The effective sample size is also affected by rows
    ## with NA values (not tested here).
    if (!use.ridge && min.subproblem < rank.required) {
        stop(gettextf("too few samples in subproblems: %.0f required, %.0f available",
                      rank.required, min.subproblem, domain = "R-sisal"),
             domain = NA)
    }

    ALL.INCR <- 100
    NAME_COEFS <- "b"
    NAME_INFO_RANK <- c("gcv.info", "rank")
    ICEPT_NAME <- "(Intercept)"

    kM <- Mtimes * kfold
    ySquared <- y2^2
    this.name <- vars.to.name(selected)
    names.old <- this.name
    n.evaluated <- rep.int(0, D)
    vertices <- character(ALL.INCR)
    vertices.logical <- matrix(FALSE, ALL.INCR, d)
    vertices[1] <- this.name
    vertices.logical[1, ] <- as.vector(selected)
    n.vertices <- 1
    n.vertex.inputs <- vector(mode = typeof(d), length = ALL.INCR)
    n.vertex.inputs[1] <- d
    E.v.all <- numeric(ALL.INCR)
    E.tr.all <- numeric(ALL.INCR)
    s.tr.all <- numeric(ALL.INCR)
    n.samples.all <- numeric(ALL.INCR)
    rank.def.all <- numeric(ALL.INCR)
    n.NA.models.all <- numeric(ALL.INCR)
    E.v.level.rank <- numeric(ALL.INCR)
    min.branches <- numeric(ALL.INCR)
    min.branches[1] <- 1
    edges <- vector(mode = "list", length = ALL.INCR)
    edges[[1]] <- list(edges=numeric(0))
    if (verbose > 0 && verbose < 2) {
        pb <- txtProgressBar(min = 0, max = D, initial = 0)
    }
    pairwise.points <- matrix(0, d, d)
    max.n.new <- max.width * hbranches
    names.new <- character(max.n.new)
    child.id <- matrix(NA_real_, max.width, hbranches)
    chosen.sets <- 1
    n.new <- 1

    ## .Random.seed is stored in user's workspace.
    ## If .Random.seed does not exist, we create one.
    ## If the function is interrupted, .Random.seed will be restored on.exit.
    if (!exists(".Random.seed", where=1)) {
        set.seed(NULL)
        on.exit(rm(".Random.seed", pos=1))
        saved.seed <- get(".Random.seed", pos=1)
    } else {
        saved.seed0 <- get(".Random.seed", pos=1)
        on.exit(assign(".Random.seed", saved.seed0, pos=1))
        seedOK <- TRUE
        tryCatch(sample(2), warning = function(...) seedOK <<- FALSE)
        if (seedOK) {
            saved.seed <- saved.seed0
            assign(".Random.seed", saved.seed, pos=1)
        } else {
            set.seed(NULL)
            saved.seed <- get(".Random.seed", pos=1)
        }
    }

    ## Check if subsets of columns of X are constant zero.  Elements
    ## of 'zero.col' are NA, FALSE or TRUE.  NA means that the results
    ## (i.e. whether the column is constant zero or not) depends on
    ## the set of variables included in the model: all samples with an
    ## NA value in any of the variables are dropped from the model,
    ## and this depends on the choice of variables.
    zero.col <- array(FALSE, c(d, Mtimes, kfold))
    for (count in seq_len(Mtimes)) {
        shuffled <- sample.int(n)
        for (slice in seq_len(kfold)) {
            idx.1 <- part.start[slice]
            idx.2 <- part.end[slice]
            idx.train <- shuffled[c(seq_len(idx.1 - 1),
                                    seq(from = idx.2 + 1,
                                        by = 1,
                                        length.out = n - idx.2))]
            X.train <- X2[idx.train, , drop=FALSE]
            zeroRange.all <- which(apply(X.train, 2, zeroRange))
            zero.col.all <- rep.int(FALSE, d)
            zero.col.all[zeroRange.all] <-
                vapply(.colMeans(X.train[, zeroRange.all, drop = FALSE],
                                 length(idx.train), length(zeroRange.all),
                                 na.rm=TRUE),
                       isZero, TRUE)
            idx.cleantrain <- idx.clean[idx.train]
            if (all(idx.cleantrain)) {
                zero.col[, count, slice] <- zero.col.all
            } else {
                zero.col.tmp <- rep.int(NA, d)
                zero.col.tmp[zero.col.all] <- TRUE
                X.noNA <- X.train[idx.cleantrain, , drop=FALSE]
                nX.noNA <- nrow(X.noNA)
                if (nX.noNA >= 2) {
                    zeroRange.noNA <- which(apply(X.noNA, 2, zeroRange))
                    zero.col.noNA <- rep.int(FALSE, d)
                    zero.col.noNA[zeroRange.noNA] <-
                        vapply(.colMeans(X.noNA[, zeroRange.noNA, drop=FALSE],
                                         length(idx.cleantrain),
                                         length(zeroRange.noNA)),
                               isZero, TRUE)
                    zero.col.tmp[!zero.col.noNA] <- FALSE
                } else if (nX.noNA == 1) {
                    zero.col.noNA <- vapply(drop(X.noNA), isZero, TRUE)
                    zero.col.tmp[!zero.col.noNA] <- FALSE
                }
                zero.col[, count, slice] <- zero.col.tmp
            }
        }
    }
    if (verbose >= 4) {
        zero.col.table <-
            apply(zero.col, 1,
                  function(x) {
                      table(factor(x, levels = as.character(c(FALSE, TRUE))),
                            useNA="always")
                  })
        colnames(zero.col.table) <- X.names
        splitcat(paste0(gettext(c("Is a column all-zero?",
                                  "Summary of all training sets.",
                                  "TRUE is yes, FALSE is no, NA is maybe:"),
                                domain = "R-sisal"),
                        collapse = " "))
        print(zero.col.table)
    }

    if (add.constant) {
        X2 <- cbind(1, X2)
    }

    ## Record whether branching brought any improvement compared to
    ## hbranches == 1, considering minimum E.v for each number of
    ## input variables.
    if (hbranches > 1) {
        ## branching.useful[k]
        ## is equivalent to
        ## E.v.best[k+1] < E.v.nobranch[k+1]
        branching.useful <- rep.int(NA, d)
        E.v.nobranch <- rep.int(NA_real_, D)
        E.tr.nobranch <- rep.int(NA_real_, D)
        s.tr.nobranch <- rep.int(NA_real_, D)
        names.nobranch <- rep.int(NA_character_, D)
    } else {
        branching.useful <- NULL
        E.v.nobranch <- NULL
        E.tr.nobranch <- NULL
        s.tr.nobranch <- NULL
    }

    while (n.selected >= 0) {
        oneSelected <- n.selected == 1
        N.selected <- n.selected + 1    # + 1 to account for the empty model
        n.coefs <- n.selected + add.constant
        n.old <- length(names.old)
        E.v <- numeric(n.old)
        E.tr <- numeric(n.old)
        s.tr <- numeric(n.old)
        n.samples <- numeric(n.old)
        rank.def <- numeric(n.old)
        n.NA.models <- numeric(n.old)
        real.branches <- min(n.selected, hbranches)
        n.new.previous <- n.new
        n.new <- 0
        n.vertices.0 <- n.vertices
        anyCoefs <- n.coefs > 0
        for (k in seq_len(n.old)) {
            this.name <- names.old[k]
            selected <- name.to.vars(this.name, d)
            constant.selected <- constant.X[selected]
            n.constants.selected <- sum(constant.selected)
            if (use.ridge && n.selected >= 1 &&
                n.constants.selected < n.selected) {
                magic.sp <- sp
                diagS <- rep.int(1, n.selected)
                ## Don't penalize constant variables
                diagS[constant.selected] <- 0
                magic.S <- list(diag(diagS))
                S.rank <- n.selected - n.constants.selected
                magic.off <- 1 + add.constant
                use.ridge2 <- TRUE
            } else {
                magic.sp <- NULL
                magic.S <- list()
                S.rank <- NULL
                magic.off <- numeric(0)
                use.ridge2 <- FALSE
            }
            if (add.constant) {
                rank.required <- n.selected + 1
            } else if (any(constant.X[selected])) {
                rank.required <- n.selected + 1 - sum(constant.X[selected])
            } else {
                rank.required <- n.selected
            }
            selected.1 <- selected
            if (add.constant) {
                selected.1 <- c(TRUE, selected.1)
            }
            n.evaluated[N.selected] <- n.evaluated[N.selected] + 1
            if (verbose >= 2) {
                cat("\n", sum(n.evaluated), " (", sum(selected), "): ",
                    paste0(which(selected), collapse=" "), "\n", sep="")
            }

            goodX <- goodX0
            tmpGoodX <- rowSums(xNA[notGood0, selected, drop=FALSE]) == 0
            goodX[notGood0] <- tmpGoodX
            n.samples[k] <- nGood0 + sum(tmpGoodX)

            ## The same Mtimes * kfold training and validation sets
            ## are used with every set of selected variables.  In
            ## theory, this ensures that training error is a
            ## decreasing function of the number of variables in case
            ## of nested models.  Computing the same randomized
            ## orderings repeatedly is quite cheap.  The alternative
            ## would be to save the orderings and reuse them.

            assign(".Random.seed", saved.seed, pos=1)
            if (verbose >= 5) {
                splitcat(paste(".Random.seed[1:6] is",
                               paste0(get(".Random.seed", pos=1,
                                          mode="numeric")[1:6],
                                      collapse=" ")))
            }
            if (anyCoefs) {
                res <- matrix(0, kM, n.coefs + 3)
                slice2 <- 0
                for (count in seq_len(Mtimes)) {
                    shuffled <- sample.int(n)
                    for (slice in seq_len(kfold)) {
                        slice2 <- slice2 + 1
                        idx.1 <- part.start[slice]
                        idx.2 <- part.end[slice]
                        idx.valid <- shuffled[idx.1:idx.2]
                        idx.valid <- idx.valid[goodX[idx.valid]]
                        idx.train <- shuffled[c(seq_len(idx.1 - 1),
                                                seq(from = idx.2 + 1,
                                                    by = 1,
                                                    length.out = n - idx.2))]
                        idx.train <- idx.train[goodX[idx.train]]
                        nGood <- length(idx.train)
                        ## It's silly to use only 1 sample, but that
                        ## is the absolute technical limit for
                        ## magic(): 0 samples will produce an
                        ## error. Cases of rank deficiency (after
                        ## regularization, if any) will be recorded.
                        if (nGood > 0) {
                            X.train <- X2[idx.train, selected.1, drop=FALSE]
                            zero.col2 <- zero.col[selected, count, slice]
                            for (colNo in which(is.na(zero.col2))) {
                                thisCol <- X.train[, colNo]
                                if (nGood > 1) {
                                    meanThis <- mean(thisCol)
                                    zero.col2[colNo] <-
                                        isTRUE(all.equal(0, meanThis)) &&
                                            zeroRange(thisCol, mean.x = meanThis)
                                } else {
                                    zero.col2[colNo] <-
                                        isTRUE(all.equal(0, thisCol))
                                }
                            }
                            nZeros <- sum(zero.col2)
                            if (nZeros < n.coefs) {
                                X.valid <- X2[idx.valid, selected.1, drop=FALSE]
                                if (add.constant) {
                                    nonzero.col <- c(TRUE, !zero.col2)
                                } else {
                                    nonzero.col <- !zero.col2
                                }
                                rank.required2 <- rank.required - nZeros
                                magic.S2 <- magic.S
                                S.rank2 <- S.rank
                                magic.sp2 <- magic.sp
                                if (nZeros > 0) {
                                    ## Constant zero columns are dropped
                                    ## from the model
                                    X.valid <-
                                        X.valid[, nonzero.col, drop=FALSE]
                                    X.train <-
                                        X.train[, nonzero.col, drop=FALSE]
                                    if (use.ridge2) {
                                        S.rank2 <- S.rank - nZeros
                                        if (S.rank2 > 0) {
                                            diagS2 <- diagS[!zero.col2]
                                            magic.S2 <- list(diag(diagS2))
                                        } else {
                                            magic.S2 <- list()
                                            magic.sp2 <- NULL
                                            S.rank2 <- NULL
                                        }
                                    }
                                }
                                y.train <- y2[idx.train]
                                the.lm <- magic(y=y.train, X=X.train,
                                                sp=magic.sp2, S=magic.S2,
                                                rank=S.rank2,
                                                off=magic.off, ...)
                                coefsTmp <- the.lm[[NAME_COEFS]]
                                coefs <- rep.int(0, N.selected)
                                coefs[nonzero.col] <- coefsTmp
                                res[slice2, 1:n.coefs] <- coefs
                                lmPredict <- X.valid %*% coefsTmp
                                lmFit <- X.train %*% coefsTmp
                                res[slice2, n.coefs + 1] <-
                                    mean((lmPredict - y2[idx.valid])^2)
                                res[slice2, n.coefs + 2] <-
                                    mean((lmFit - y.train)^2)
                                res[slice2, n.coefs + 3] <-
                                    !isTRUE(the.lm[[NAME_INFO_RANK]] >=
                                            rank.required2)
                            } else {
                                ## All variables are constant zeros,
                                ## get zero coefficients
                                res[slice2, 1:n.coefs] <- 0
                                res[slice2, n.coefs + 1] <-
                                    mean(ySquared[idx.valid])
                                res[slice2, n.coefs + 2] <-
                                    mean(ySquared[idx.train])
                                res[slice2, n.coefs + 3] <- FALSE
                            }
                        } else {
                            ## No samples available => all variables
                            ## get NA coefficients
                            res[slice2, 1:n.coefs] <- NA_real_
                            res[slice2, n.coefs + 1] <-
                                mean(ySquared[idx.valid])
                            res[slice2, n.coefs + 2] <-
                                mean(ySquared[idx.train])
                            res[slice2, n.coefs + 3] <- TRUE
                        }
                    }
                }
                rank.def[k] <- sum(res[, n.coefs + 3])
                mse.train <- res[, n.coefs + 2]
                mse.valid <- res[, n.coefs + 1]
                b <- res[, 1:n.coefs, drop=FALSE]
            } else {
                mse.valid <- numeric(kM)
                mse.train <- numeric(kM)
                slice2 <- 0
                for (count in seq_len(Mtimes)) {
                    shuffled <- sample.int(n)
                    for (slice in seq_len(kfold)) {
                        slice2 <- slice2 + 1
                        idx.1 <- part.start[slice]
                        idx.2 <- part.end[slice]
                        idx.valid <- shuffled[idx.1:idx.2]
                        idx.train <- shuffled[c(seq_len(idx.1 - 1),
                                                seq(from = idx.2 + 1,
                                                    by = 1,
                                                    length.out = n - idx.2))]
                        mse.valid[slice2] <- mean(ySquared[idx.valid])
                        mse.train[slice2] <- mean(ySquared[idx.train])
                    }
                }
            }

            E.tr[k] <- mean(mse.train)
            s.tr[k] <- sd(mse.train)
            this.E.v <- mean(mse.valid, na.rm=TRUE)
            E.v[k] <- this.E.v
            if (verbose >= 3) {
                splitcat(gettextf("E.v is %f", this.E.v, domain = "R-sisal"))
            }
            if (n.selected > 0) {
                n.NA.models[k] <- sum(is.na(b[, 1]))
                ## There is always at least one non-NA row in b (see
                ## test above which checks for presence of at least 2
                ## complete samples)
                m.b <- apply(b, 2, median, na.rm = TRUE)
                q.b <- apply(b, 2, quantile, probs=c(q, 1-q), na.rm = TRUE)
                d.b <- q.b[2, ] - q.b[1, ]
                d.b.zero <- d.b == 0
                ratio.b <- abs(m.b) / d.b

                ## very unlikely
                if (any(d.b.zero)) {
                    ## even more unlikely
                    if (all(d.b.zero)) {
                        if (add.constant) {
                            order.ratio <- sort.list(abs(m.b[-1]))
                        } else {
                            order.ratio <- sort.list(abs(m.b))
                        }
                    } else {
                        ## Treat zero deviance as minimum non-zero deviance
                        ratio.b[d.b.zero] <-
                            abs(m.b[d.b.zero]) / min(d.b[!d.b.zero])
                        if (add.constant) {
                            mask <- c(FALSE, rep.int(TRUE, n.selected))
                            idx.notzero <- which(mask & !d.b.zero)
                            idx.zero <- which(mask & d.b.zero)
                        } else {
                            idx.notzero <- which(!d.b.zero)
                            idx.zero <- which(d.b.zero)
                        }
                        ## In case of a tie between the ratio.b values of
                        ## a non-zero deviance and a zero deviance
                        ## variable, this stable sort means that the zero
                        ## deviance variable comes later in order.ratio
                        ## (the non-zero deviance variable may be dropped
                        ## and the zero deviance variable saved)
                        order.ratio <-
                            c(idx.notzero,
                              idx.zero)[sort.list(c(ratio.b[idx.notzero],
                                                    ratio.b[idx.zero]))]
                    }
                } else if (add.constant) {
                    order.ratio <- sort.list(ratio.b[-1])
                } else {
                    order.ratio <- sort.list(ratio.b)
                }
                which.selected <- which(selected)
                if (verbose >= 4) {
                    splitcat(gettext("Variables by importance, increasing order:",
                                     domain = "R-sisal"))
                    print(which.selected[order.ratio])
                    splitcat("abs(median(b)) / spread(b):")
                    if (add.constant) {
                        print(ratio.b[order.ratio + 1])
                    } else {
                        print(ratio.b[order.ratio])
                    }
                }
                edges.this <- numeric(real.branches)
                idx.vert.k <- n.vertices.0 - n.new.previous + chosen.sets[k]
                for (l in seq_len(real.branches)) {
                    candidate <- selected
                    rm.idx <- which.selected[order.ratio[l]]
                    if (verbose >= 4) {
                        splitcat(gettextf("Dropping variable %.0f", rm.idx,
                                          domain = "R-sisal"))
                    }
                    candidate[rm.idx] <- FALSE
                    ## Remaining variables are "better", except those that
                    ## were actually dropped earlier, i.e. in another branch
                    if (n.selected >= 2) {
                        better <- candidate
                        better[which.selected[order.ratio[seq_len(l - 1)]]] <-
                            FALSE
                        pairwise.points[better, rm.idx] <-
                            pairwise.points[better, rm.idx] + 1
                    }
                    candidate.name <- vars.to.name(candidate, oneSelected)
                    candidate.id <-
                        which(names.new[seq_len(n.new)] == candidate.name)[1]
                    if (is.na(candidate.id)) {
                        n.new <- n.new + 1
                        candidate.id <- n.new
                        names.new[candidate.id] <- candidate.name
                        if(verbose >= 3) {
                            splitcat(gettextf("New candidate %s",
                                              candidate.name,
                                              domain = "R-sisal"))
                        }
                        ## The vectors and the matrix grow in chunks
                        ## of ALL.INCR elements or rows
                        length.now <- length(vertices)
                        if (length.now == n.vertices) {
                            length.after <- length.now + ALL.INCR
                            length(vertices) <- length.after
                            vertices.logical <-
                                rbind(vertices.logical,
                                      matrix(FALSE, ALL.INCR, d))
                            length(n.vertex.inputs) <- length.after
                            length(edges) <- length.after
                            length(min.branches) <- length.after
                        }
                        n.vertices <- n.vertices + 1
                        vertices[n.vertices] <- candidate.name
                        vertices.logical[n.vertices, ] <- as.vector(candidate)
                        n.vertex.inputs[n.vertices] <- n.selected - 1
                        edges[[n.vertices]] <- list(edges=numeric(0))
                        min.branches[n.vertices] <-
                            max(l, min.branches[idx.vert.k])
                    } else if (verbose >= 3) {
                        splitcat(gettextf("Duplicate candidate %s",
                                          candidate.name, domain = "R-sisal"))
                        idx.vert.candidate <- n.vertices.0 + candidate.id
                        min.branches[idx.vert.candidate] <-
                            min(min.branches[idx.vert.candidate],
                                max(l, min.branches[idx.vert.k]))
                    }
                    child.id[k, l] <- candidate.id
                    edges.this[l] <- n.vertices.0 + candidate.id
                }
                edges[[idx.vert.k]][[1]] <- edges.this
            }
        }
        ## Stable sort
        level.order.E.v <- sort.list(E.v)
        E.v.rank.this <- numeric(n.old)
        E.v.rank.this[level.order.E.v] <- 1:n.old
        ## In case of ties (unlikely), nodes evaluated earlier are
        ## preferred.  If pruning.keep.best is TRUE, the first node
        ## evaluated on each round is the no-branching node.
        ## Therefore, if the node chosen as the best one is different
        ## than the no-branching node, the error of the best node is
        ## surely smaller than that of the no-branching node. NOTE:
        ## This holds when the nodes have the same number of
        ## variables.
        which.best <- level.order.E.v[1]
        idx.append <- n.vertices.0 - n.new.previous + chosen.sets
        if (hbranches > 1) {
            if (pruning.keep.best) {
                idx.nobranch <- 1
            } else {
                idx.nobranch <- match(1, min.branches[idx.append])
            }
            if (!is.na(idx.nobranch)) {
                if (which.best != idx.nobranch) {
                    ## Branch not taken when n.selected == 0
                    branching.useful[n.selected] <- TRUE
                } else if (n.selected > 0) {
                    branching.useful[n.selected] <- FALSE
                }
                E.v.nobranch[N.selected] <- E.v[idx.nobranch]
                s.tr.nobranch[N.selected] <- s.tr[idx.nobranch]
                E.tr.nobranch[N.selected] <- E.tr[idx.nobranch]
                names.nobranch[N.selected] <- names.old[idx.nobranch]
            }
        }
        E.v.best[N.selected] <- E.v[which.best]
        s.tr.best[N.selected] <- s.tr[which.best]
        E.tr.best[N.selected] <- E.tr[which.best]
        names.best[N.selected] <- names.old[which.best]
        length.now <- length(E.v.all)
        length.diff <- n.vertices.0 - length.now
        ## Vectors grow in chunks that are multiple of ALL.INCR in size
        if (length.diff > 0) {
            length.diff <- ceiling(length.diff / ALL.INCR) * ALL.INCR
            length.after <- length.now + length.diff
            length(E.v.all) <- length.after
            length(E.tr.all) <- length.after
            length(s.tr.all) <- length.after
            length(rank.def.all) <- length.after
            length(n.NA.models.all) <- length.after
            length(E.v.level.rank) <- length.after
        }
        E.v.all[idx.append] <- E.v
        E.tr.all[idx.append] <- E.tr
        s.tr.all[idx.append] <- s.tr
        n.samples.all[idx.append] <- n.samples
        E.v.level.rank[idx.append] <- E.v.rank.this
        rank.def.all[idx.append] <- rank.def
        n.NA.models.all[idx.append] <- n.NA.models
        if (verbose >= 3) {
            cat("\n")
            if (n.new <= max.width) {
                splitcat(gettextf("Not pruning (%d <= %d)",
                                  n.new, max.width, domain = "R-sisal"))
            } else {
                splitcat(gettextf("Pruning (%d > %d)...",
                                  n.new, max.width, domain = "R-sisal"))
            }
        }
        if (n.new <= max.width) {
            chosen.sets <- seq_len(n.new)
        } else if (prune.rr || prune.greedy) { # round robin or greedy
            if (pruning.reverse) {
                order.E.v <- sort.list(E.v, decreasing = TRUE)
                order.branches <- real.branches:1
            } else {
                order.E.v <- sort.list(E.v)
                order.branches <- 1:real.branches
            }
            if (prune.rr) {
                chosen.sets <- as.numeric(child.id[order.E.v, order.branches])
            } else {
                chosen.sets <- as.numeric(t(child.id[order.E.v,
                                                     order.branches]))
            }
            chosen.idx <- which(!duplicated(chosen.sets))
            idx.1 <- which(chosen.sets[chosen.idx] == 1)
            if (idx.1 <= max.width && idx.1 > 1) {
                chosen.idx <- chosen.idx[c(idx.1, seq_len(idx.1-1),
                                           seq(from = idx.1 + 1, by = 1,
                                               length.out=max.width - idx.1))]
            } else if (idx.1 > max.width && pruning.keep.best) {
                chosen.idx <- chosen.idx[c(idx.1, seq_len(max.width - 1))]
            } else {
                chosen.idx <- chosen.idx[seq_len(max.width)]
            }
            chosen.sets <- chosen.sets[chosen.idx]
            if (verbose >= 3) {
                if (prune.rr) {
                    parent.id <- rep.int(order.E.v, real.branches)[chosen.idx]
                    child.rank <- rep(order.branches, each=n.old)[chosen.idx]
                } else {
                    parent.id <- rep(order.E.v, each=real.branches)[chosen.idx]
                    child.rank <- rep.int(order.branches, n.old)[chosen.idx]
                }
                print(data.frame(Parent=parent.id, Rank=child.rank,
                                 ID=chosen.sets))
            }
        } else if (prune.edges) { # random edges
            ## Get fresh random numbers
            if (exists("saved.seed2", inherits = FALSE)) {
                assign(".Random.seed", saved.seed2, pos=1)
            }
            if (verbose >= 5) {
                splitcat(paste(".Random.seed[1:6] is",
                               paste0(get(".Random.seed", pos=1,
                                          mode="numeric")[1:6],
                                      collapse=" ")))
            }
            child.counts <-
                table(as.numeric(child.id[1:n.old, 1:real.branches]))
            if (pruning.keep.best) {
                child.counts <- child.counts[-which(names(child.counts) == "1")]
            }
            sample.from <- as.numeric(names(child.counts))
            if (pruning.reverse) {
                sample.prob <- 1 / as.numeric(child.counts)
            } else {
                sample.prob <- as.numeric(child.counts)
            }
            if (pruning.keep.best) {
                chosen.sets <-
                    c(1, sample(x = sample.from, size = max.width - 1,
                                prob = sample.prob, replace = FALSE))
            } else {
                chosen.sets <- sample(x = sample.from, size = max.width,
                                      prob = sample.prob, replace = FALSE)
            }
            if (verbose >= 4) {
                splitcat(gettext("Sampling from:", domain = "R-sisal"))
                print(sample.from)
                splitcat(gettext("with unnormalized probabilities:",
                                 domain = "R-sisal"))
                print(sample.prob)
            }
            if (verbose >= 3) {
                cat(gettext("IDs:\n", domain = "R-sisal"))
                print(chosen.sets)
            }
            saved.seed2 <- get(".Random.seed", pos=1, mode="numeric")
        } else { # "random nodes"
            ## Get fresh random numbers
            if (exists("saved.seed2", inherits = FALSE)) {
                assign(".Random.seed", saved.seed2, pos=1)
            }
            if (verbose >= 5) {
                splitcat(paste(".Random.seed[1:6] is",
                               paste0(get(".Random.seed", pos=1,
                                          mode="numeric")[1:6],
                                      collapse=" ")))
            }
            if (pruning.keep.best) {
                chosen.sets <-
                    c(1, sample(x = seq(from = 2, to = n.new),
                                size = max.width - 1, replace=FALSE))
            } else {
                chosen.sets <- sample.int(n = n.new,
                                          size = max.width, replace=FALSE)
            }
            if (verbose >= 3) {
                cat(gettext("IDs:\n", domain = "R-sisal"))
                print(chosen.sets)
            }
            saved.seed2 <- get(".Random.seed", pos=1, mode="numeric")
        }
        names.old <- names.new[chosen.sets]
        if (exists("pb", inherits = FALSE)) {
            setTxtProgressBar(pb, value = D - n.selected)
        }
        n.selected <- n.selected - 1
    }
    names(E.v.best) <- names.best
    names(s.tr.best) <- names.best
    names(E.tr.best) <- names.best
    idx.L.v <- which.min(E.v.best)
    n.L.v <- idx.L.v - 1
    L.v <- name.to.vars.idx(names.best[idx.L.v])
    names(L.v) <- X.names[L.v]
    idx.L.f <- min(which(E.v.best <= E.v.best[idx.L.v] + s.tr.best[idx.L.v]))
    n.L.f <- idx.L.f - 1
    if (n.L.f < n.L.v) {
        L.f <- name.to.vars.idx(names.best[idx.L.f])
        names(L.f) <- X.names[L.f]
    } else {
        L.f <- L.v
    }
    if (hbranches > 1) {
        names(E.v.nobranch) <- names.nobranch
        names(s.tr.nobranch) <- names.nobranch
        names(E.tr.nobranch) <- names.nobranch
    }
    if (!is.null(branching.useful) && !anyNA(branching.useful)) {
        idx.L.v.nobranch <- which.min(E.v.nobranch)
        n.L.v.nobranch <- idx.L.v.nobranch - 1
        L.v.nobranch <- name.to.vars.idx(names.nobranch[idx.L.v.nobranch])
        names(L.v.nobranch) <- X.names[L.v.nobranch]
        idx.L.f.nobranch <- min(which(E.v.nobranch <=
                                      E.v.nobranch[idx.L.v.nobranch] +
                                      s.tr.nobranch[idx.L.v.nobranch]))
        n.L.f.nobranch <- idx.L.f.nobranch - 1
        if (n.L.f.nobranch < n.L.v.nobranch) {
            L.f.nobranch <- name.to.vars.idx(names.nobranch[idx.L.f.nobranch])
            names(L.f.nobranch) <- X.names[L.f.nobranch]
        } else {
            L.f.nobranch <- L.v.nobranch
        }
    } else {
        L.v.nobranch <- NULL
        L.f.nobranch <- NULL
    }
    idx.vert <- seq_len(n.vertices)
    vertices <- vertices[idx.vert]
    vertices.logical <- vertices.logical[idx.vert, , drop=FALSE]
    n.vertex.inputs <- n.vertex.inputs[idx.vert]
    edges <- edges[idx.vert]
    names(edges) <- vertices

    ## Ridge regression wrapper to magic()
    magicRidge <- function(X, y, sp, add.constant, ...) {
        if (ncol(X) == 0 && !add.constant) {
            return(NULL)
        }
        idx.clean <- rowSums(is.na(X)) == 0 & !is.na(y)
        X2 <- X[idx.clean, , drop = FALSE]
        y2 <- y[idx.clean]
        constant.X <- apply(X2, 2, zeroRange)
        d <- ncol(X)
        ## Constant variables are not penalized
        nConstants <- sum(constant.X)
        if (nConstants == d) {
            sp2 <- NULL
            S <- list()
            S.rank <- NULL
            off <- numeric(0)
        } else {
            sp2 <- sp
            diagS <- rep.int(1, d)
            diagS[constant.X] <- 0
            S <- list(diag(diagS))
            S.rank <- d - nConstants
            if (add.constant) {
                off <- 2
            } else {
                off <- 1
            }
        }
        if (add.constant) {
            X2 <- cbind(rep.int(1, nrow(X2)), X2)
        }
        magic(y = y2, X = X2,
              sp = sp2, S = S, rank = S.rank, off = off, ...)
    }

    if (add.constant) {
        lmForm <- y.complete ~ .
        iceptName <- ICEPT_NAME
    } else {
        lmForm <- y.complete ~ . - 1
        iceptName <- character(0)
    }
    lm.full <- lm(lmForm, data=as.data.frame(X.complete),
                  na.action=na.exclude)
    magic.full <- magicRidge(X = X.complete, y = y.complete,
                             add.constant = add.constant, sp = sp, ...)
    names(magic.full[[NAME_COEFS]]) <- c(iceptName, X.names)
    if (n.L.v < d) {
        X.L.v <- X.complete[, L.v, drop=FALSE]
        if (n.L.v > 0) {
            lm.L.v <- lm(lmForm, data=as.data.frame(X.L.v),
                         na.action=na.exclude)
        } else if (add.constant) {
            lm.L.v <- lm(y.complete ~ 1, na.action=na.exclude)
        } else {
            lm.L.v <- NULL
        }
        magic.L.v <- magicRidge(X = X.L.v, y = y.complete,
                                add.constant = add.constant, sp = sp, ...)
        if (!is.null(magic.L.v)) {
            names(magic.L.v[[NAME_COEFS]]) <- c(iceptName, X.names[L.v])
        }
    } else {
        lm.L.v <- lm.full
        magic.L.v <- magic.full
    }
    if (n.L.f < n.L.v) {
        X.L.f <- X.complete[, L.f, drop=FALSE]
        if (n.L.f > 0) {
            lm.L.f <- lm(lmForm, data=as.data.frame(X.L.f),
                         na.action=na.exclude)
        } else if (add.constant) {
            lm.L.f <- lm(y.complete ~ 1, na.action=na.exclude)
        } else {
            lm.L.f <- NULL
        }
        magic.L.f <- magicRidge(X = X.L.f, y = y.complete,
                                add.constant = add.constant, sp = sp, ...)
        if (!is.null(magic.L.f)) {
            names(magic.L.f[[NAME_COEFS]]) <- c(iceptName, X.names[L.f])
        }
    } else {
        lm.L.f <- lm.L.v
        magic.L.f <- magic.L.v
    }

    E.v.level.rank <- E.v.level.rank[idx.vert]
    rank.def.all <- rank.def.all[idx.vert]
    n.NA.models.all <- n.NA.models.all[idx.vert]
    vertex.data <- data.frame(E.tr = E.tr.all[idx.vert],
                              s.tr = s.tr.all[idx.vert],
                              E.v = E.v.all[idx.vert],
                              n.samples = n.samples.all[idx.vert],
                              E.v.level.rank = E.v.level.rank,
                              n.rank.deficient = rank.def.all,
                              n.NA.models = n.NA.models.all,
                              n.inputs = n.vertex.inputs,
                              min.branches = min.branches[idx.vert],
                              row.names = vertices)
    colnames(vertices.logical) <- X.names
    rank.def.nodes <- sum(rank.def.all > 0)
    if (rank.def.nodes > 0) {
        tmp <-
            substitute(sprintf(ngettext(rank.def.nodes,
                                        "%.0f set of inputs with rank deficiency problems",
                                        "%.0f sets of inputs with rank deficiency problems",
                                        domain = "R-sisal"),
                               rank.def.nodes))
        n.warn <- n.warn + 1
        if (n.warn <= max.warn) {
            allWarn[[n.warn]] <- tmp
        }
        warning(eval(tmp), domain = NA)
        n.NA.tmp <- as.vector(na.omit(n.NA.models.all))
        NA.nodes <- sum(n.NA.tmp > 0)
        if (NA.nodes > 0) {
            tmp <-
                substitute(sprintf(ngettext(NA.nodes,
                                            "%.0f set of inputs with NA model coefficients",
                                            "%.0f sets of inputs with NA model coefficients",
                                            domain = "R-sisal"),
                                   NA.nodes))
            n.warn <- n.warn + 1
            if (n.warn <= max.warn) {
                allWarn[[n.warn]] <- tmp
            }
            warning(eval(tmp), domain = NA)
            if (verbose >= 4) {
                NA.percentage <- n.NA.tmp / kM * 100
                splitcat(paste0(gettext("Percentage of NA models per node",
                                        "quantiles of values > 0:",
                                        domain = "R-sisal"),
                                collapse = "; "))
                print(quantile(NA.percentage[NA.percentage > 0]))
            }
        }
    }
    ## Copeland's pairwise aggregation method.  In this (classic?)
    ## formulation, the score of a candidate is the number of pairwise
    ## victories minus the number of pairwise defeats.  So, win = +1,
    ## tie = 0, loss = -1.  In another formulation, win = +1, tie =
    ## +1/2, loss = 0.  We see these are equivalent.
    pairwise.wins <- pairwise.points > t(pairwise.points)
    dimnames(pairwise.points) <-
        list("times removed later than" = X.names,
             "this other input" = X.names)
    dimnames(pairwise.wins) <-
        list("is more important than" = X.names,
             "this other input" = X.names)
    ## wins - losses (victories - defeats)
    pairwise.preferences <- rowSums(pairwise.wins) - colSums(pairwise.wins)
    pairwise.rank <- rank(-pairwise.preferences, ties.method = "min")

    ## Examine the best nodes for each number of inputs
    rank1.vd <- vertex.data[which(E.v.level.rank == 1), , drop=FALSE]
    vd.order <- sort.list(rank1.vd[["n.inputs"]], decreasing=TRUE)
    rank1.vd <- rank1.vd[vd.order, , drop=FALSE]
    rank1.names <- row.names(rank1.vd)
    rank1.sets <- lapply(rank1.names, name.to.vars.idx)
    path.length <- 1
    ## nested.path is the removal order of variables on a path
    ## containing the incrementally smaller best models with a given
    ## number of variables, if the models are nested.  If the models
    ## are not nested, the value of nested.path is invalid.  Even in
    ## the case that the models are nested, the path may not actually
    ## exist in the search graph, but that does not matter for the
    ## purpose of ranking.
    nested.path <- numeric(d)
    for (k in seq_len(d)) {
        this.diff <- setdiff(rank1.sets[[k]], rank1.sets[[k + 1]])
        if (length(this.diff) == 1) {
            n.paths <- length(path.length)
            path.length[n.paths] <- path.length[n.paths] + 1
            nested.path[k] <- this.diff
        } else {
            path.length <- c(path.length, 1)
        }
    }
    if (length(path.length) == 1) {
        ## Add the remaining input variable to the path
        nested.rank <- numeric(d)
        nested.rank[nested.path] <- d : 1
        names(nested.rank) <- X.names
        names(nested.path) <- X.names[nested.path]
    } else {
        nested.path <- NULL
        nested.rank <- NULL
    }
    if (n.warn < max.warn) {
        allWarn <- allWarn[seq_len(n.warn)]
    }

    res <- list(L.f=L.f, L.v=L.v, E.tr=E.tr.best, s.tr=s.tr.best, E.v=E.v.best,
                L.f.nobranch=L.f.nobranch, L.v.nobranch=L.v.nobranch,
                E.tr.nobranch=E.tr.nobranch, s.tr.nobranch=s.tr.nobranch,
                E.v.nobranch=E.v.nobranch,
                n.evaluated=n.evaluated, edges=edges, vertices=vertices,
                vertices.logical=vertices.logical,
                vertex.data=vertex.data, var.names=orig.names,
                n=n, d=d, n.missing=n.missing, n.clean=n.clean,
                lm.L.f=lm.L.f, lm.L.v=lm.L.v, lm.full=lm.full,
                magic.L.f=magic.L.f, magic.L.v=magic.L.v, magic.full=magic.full,
                mean.y=mean.y, sd.y=sd.y, zeroRange.y=zeroRange.y,
                mean.X=mean.X, sd.X=sd.X, zeroRange.X=zeroRange.X,
                constant.X=constant.X, params=params,
                pairwise.points=pairwise.points,
                pairwise.wins=pairwise.wins,
                pairwise.preferences=pairwise.preferences,
                pairwise.rank=pairwise.rank,
                path.length=path.length, nested.path=nested.path,
                nested.rank=nested.rank, branching.useful=branching.useful,
                warnings=allWarn, n.warn=n.warn)
    class(res) <- c("sisal")
    on.exit() # if the function returns uninterrupted, cancel cleanup measures
    res
}

print.sisal <- function(x, max.warn=10, ...) {
    if (!inherits(x, "sisal")) {
        stop('use only with "sisal" objects')
    }
    stopifnot(is.numeric(max.warn), length(max.warn) == 1,
              round(max.warn) == max.warn, max.warn >= 0)
    params <- x[["params"]]
    hbranches <- params[["hbranches"]]
    L.v <- x[["L.v"]]
    L.f <- x[["L.f"]]
    n.L.v <- length(L.v)
    n.L.f <- length(L.f)
    sisaltype <- paste(ifelse(hbranches > 1,
                              gettext("branching", domain="R-sisal"),
                              gettext("no branching", domain="R-sisal")),
                       ifelse(params[["use.ridge"]],
                              gettext("ridge regression", domain="R-sisal"),
                              gettext("OLS regression", domain="R-sisal")),
                       sep=", ")
    splitcat(gettextf("sisal (%s) results", sisaltype, domain="R-sisal"))

    cat("\n")
    cat(gettext("Parameters:\n", domain="R-sisal"))
    splitcat(paste(names(params), params, sep="=", collapse=", "))

    cat("\n")
    d <- x[["d"]]
    n <- x[["n"]]
    splitcat(gettext("Data dimensions:", domain="R-sisal"),
             gettextf("inputs: %.0f, samples: %.0f", d, n,
                      domain="R-sisal"))

    cat("\n")
    var.names <- x[["var.names"]]
    if (is.null(var.names)) {
        splitcat(gettext("Input names are absent", domain="R-sisal"))
        var.names.new <- paste0("V", seq_len(d))
    } else {
        splitcat(gettext("Input names are present:", domain="R-sisal"))
        print(var.names, ...)
    }

    cat("\n")
    splitcat(gettextf("Selected inputs (%.0f), model with smallest validation error (L.v):",
                      n.L.v, domain="R-sisal"))
    print(L.v, ...)

    cat("\n")
    splitcat(gettextf("Selected inputs (%.0f), least complex model with error inside threshold (L.f):",
                      n.L.f, domain="R-sisal"))
    print(L.f, ...)


    cat("\n")
    if (n.L.f > 0) {
        if (all(L.f %in% L.v)) {
            splitcat(gettext("L.f is a subset of L.v", domain="R-sisal"))
        } else {
            splitcat(gettext("L.f is not a subset of L.v", domain="R-sisal"))
        }
    }

    cat("\n")

    if (length(x[["path.length"]]) == 1 && hbranches > 1) {
        splitcat(gettext('The best input sets for each number of inputs are nested.',
                         domain="R-sisal"))
        rankingNested <- x[["nested.rank"]]
        rankingPair <- x[["pairwise.rank"]]
        if (any(rankingNested != rankingPair)) {
            splitcat(gettext("Ranking based on the removal order of inputs defining the nested input sets:",
                             domain="R-sisal"))
            if (is.null(var.names)) {
                names(rankingNested) <- var.names.new
            }
            print(rankingNested, ...)
            splitcat(gettext("Ranking based on pairwise comparison of inputs (Copeland's method):",
                             domain="R-sisal"))
        } else {
            splitcat(gettext("Ranking of inputs based on",
                             "* removal order of inputs",
                             "* pairwise comparison of inputs (Copeland's method)",
                             "Both methods agree. The ranking is:",
                             domain="R-sisal"))
        }
        if (is.null(var.names)) {
            names(rankingPair) <- var.names.new
        }
        print(rankingPair, ...)
    } else {
        ranking <- x[["pairwise.rank"]]
        if (is.null(var.names)) {
            names(ranking) <- var.names.new
        }
        if (hbranches > 1) {
            splitcat(gettext('The best input sets for each number of inputs are not nested.',
                             domain="R-sisal"))
            splitcat(gettext("Ranking based on pairwise comparison of inputs (Copeland's method):",
                             domain="R-sisal"))
        } else {
            splitcat(gettext("Removal order of input variables (ends with remaining input):",
                             domain="R-sisal"))
            print(x[["nested.path"]], ...)
            splitcat(gettext("Ranking based on the removal order:",
                             domain="R-sisal"))
        }
        print(ranking, ...)
    }

    branching.useful <- x[["branching.useful"]]
    if (!is.null(branching.useful)) {
        cat("\n")
        branchingNA <- is.na(branching.useful[seq_len(d - 1)])
        any_NA <- any(branchingNA)
        if (d == 1) {
            splitcat(gettext("The search did not branch because there was only one variable.",
                             domain="R-sisal"))
        } else if (branchingNA[d - 1]) {
            splitcat(gettext("Due to pruning of the search graph, all the \"no branching\" results are missing.",
                             domain="R-sisal"))
        } else {
            if (any_NA) {
                firstNA <- which(branchingNA)
                firstNA <- firstNA[length(firstNA)]
                splitcat(sprintf(ngettext(firstNA,
                                          "Due to pruning of the search graph, the \"no branching\" solution for %.0f variable is not available.",
                                          "Due to pruning of the search graph, the \"no branching\" solutions for %.0f or fewer variables are not available.",
                                          domain = "R-sisal"),
                                 firstNA))
                splitcat(gettext("The alternative \"no branching\" L.v and L.f sets could not be computed.",
                                 domain="R-sisal"))
            } else {
                L.v.nobranch <- x[["L.v.nobranch"]]
                L.f.nobranch <- x[["L.f.nobranch"]]
                n.L.v.nobranch <- length(L.v.nobranch)
                n.L.f.nobranch <- length(L.f.nobranch)
                E.v.best <- x[["E.v"]]
                E.v.nobranch <- x[["E.v.nobranch"]]
                equal.L.v <- FALSE
                nested.L.v <- NA
                if (n.L.v < n.L.v.nobranch) {
                    splitcat(gettextf("Branching reduced the size of the %s set.",
                                      "L.v", domain="R-sisal"))
                    nested.L.v <- all(L.v %in% L.v.nobranch)
                } else if (n.L.v > n.L.v.nobranch) {
                    splitcat(gettextf("Branching increased the size of the %s set.",
                                      "L.v", domain="R-sisal"))
                    nested.L.v <- all(L.v.nobranch %in% L.v)
                } else if (setequal(L.v, L.v.nobranch)) {
                    splitcat(gettextf("Branching did not affect the %s set.",
                                      "L.v", domain="R-sisal"))
                    equal.L.v <- TRUE
                } else {
                    splitcat(gettextf("Branching changed the composition but not the size of the %s set.",
                                      "L.v", domain="R-sisal"))
                }
                if (identical(nested.L.v, TRUE)) {
                    splitcat(gettextf("The branching and non-branching %s sets are nested.",
                                      "L.v", domain="R-sisal"))
                } else if (identical(nested.L.v, FALSE)) {
                    splitcat(gettextf("The branching and non-branching %s sets are not nested.",
                                      "L.v", domain="R-sisal"))
                }
                if (!equal.L.v) {
                    ELv.best <- E.v.best[n.L.v + 1]
                    ELv.nobranch <- E.v.nobranch[n.L.v.nobranch + 1]
                    if (ELv.best < ELv.nobranch ) {
                        splitcat(gettext("Validation error in the \"L.v, branching\" set is smaller than that in the \"L.v, no branching\" set.",
                                         domain="R-sisal"))
                    } else {
                        ## Should be a rare case
                        splitcat(gettextf("Validation error is the same in both %s sets (with or without branching)",
                                          "L.v", domain="R-sisal"))
                    }
                    splitcat(gettextf("Selected inputs (%.0f), model with smallest validation error, no branching (L.v.nobranch):",
                                      n.L.v.nobranch,, domain="R-sisal"))
                    print(L.v.nobranch, ...)
                }
                equal.L.f <- FALSE
                nested.L.f <- NA
                if (n.L.f < n.L.f.nobranch) {
                    splitcat(gettextf("Branching reduced the size of the %s set.",
                                      "L.f", domain="R-sisal"))
                    nested.L.f <- all(L.f %in% L.f.nobranch)
                } else if (n.L.f > n.L.f.nobranch) {
                    splitcat(gettextf("Branching increased the size of the %s set.",
                                      "L.f", domain="R-sisal"))
                    nested.L.f <- all(L.f.nobranch %in% L.f)
                } else if (setequal(L.f, L.f.nobranch)) {
                    splitcat(gettextf("Branching did not affect the %s set.",
                                      "L.f", domain="R-sisal"))
                    equal.L.f <- TRUE
                } else {
                    splitcat(gettextf("Branching changed the composition but not the size of the %s set.",
                                      "L.f", domain="R-sisal"))
                }
                if (identical(nested.L.f, TRUE)) {
                    splitcat(gettextf("The branching and non-branching %s sets are nested.",
                                      "L.f", domain="R-sisal"))
                } else if (identical(nested.L.f, FALSE)) {
                    splitcat(gettextf("The branching and non-branching %s sets are not nested.",
                                      "L.f", domain="R-sisal"))
                }
                if (!equal.L.f) {
                    ELf.best <- E.v.best[n.L.f + 1]
                    ELf.nobranch <- E.v.nobranch[n.L.f.nobranch + 1]
                    if (ELf.best < ELf.nobranch ) {
                        splitcat(gettext("Validation error in the \"L.f, branching\" set is smaller than that in the \"L.f, no branching\" set.",
                                         domain="R-sisal"))
                    } else if (ELf.best > ELf.nobranch) {
                        splitcat(gettext("Validation error in the \"L.f, branching\" set is larger than that in the \"L.f, no branching\" set.",
                                         domain="R-sisal"))
                    } else {
                        ## Should be a rare case
                        splitcat(gettextf("Validation error is the same in both %s sets (with or without branching)",
                                          "L.f", domain="R-sisal"))
                    }
                    splitcat(gettextf("Selected inputs (%.0f), least complex model with error inside threshold, no branching (L.f.nobranch):",
                                      n.L.f.nobranch,, domain="R-sisal"))
                    print(L.f.nobranch, ...)
                }
            }

        }
        if (any(branching.useful)) {
            splitcat(gettext("Branching reduced validation error when the number of variables was one of the following:", domain="R-sisal"))
            nvars.useful <- which(branching.useful)
            print(nvars.useful)
        } else if (any_NA) {
            splitcat(gettext("For the known part, branching did not reduce validation error.", domain="R-sisal"))
        } else {
            splitcat(gettext("Branching did not reduce validation error.",
                             domain="R-sisal"))
        }
    }

    the.warnings <- x[["warnings"]]
    ## Safeguard against arbitrary code execution
    safe.warnings <- which(vapply(the.warnings, safeToEval, TRUE))
    n.show <- min(length(safe.warnings), max.warn)
    n.warn <- x[["n.warn"]]
    if (n.warn > 0 && n.show == 0) {
        cat("\n")
        splitcat(sprintf(ngettext(n.warn,
                                  "%.0f warning from sisal(), not shown",
                                  "%.0f warnings from sisal(), not shown",
                                  domain="R-sisal"),
                         n.warn))
    } else if (n.warn > 0 && n.show < n.warn) {
        cat("\n")
        splitcat(sprintf(ngettext(n.warn,
                                  "%.0f warning from sisal(), %.0f shown:",
                                  "%.0f warnings from sisal(), %.0f shown:",
                                  domain="R-sisal"),
                         n.warn, n.show),
                 lapply(the.warnings[safe.warnings[seq_len(n.show)]], eval))
    } else if (n.warn > 0) {
        cat("\n")
        splitcat(sprintf(ngettext(n.warn,
                                  "%.0f warning from sisal():",
                                  "%.0f warnings from sisal():",
                                  domain="R-sisal"),
                         n.warn),
                 lapply(the.warnings[safe.warnings], eval))
    }
    invisible(x)
}

summary.sisal <- function(object, ...) {
    if(!inherits(object, "sisal")) {
        stop('use only with "sisal" objects')
    }
    d <- object[["d"]]
    D <- d + 1
    L.v.flag <- rep.int(FALSE, D)
    L.v.idx <- length(object[["L.v"]]) + 1
    L.v.flag[L.v.idx] <- TRUE
    L.f.flag <- rep.int(FALSE, D)
    L.f.idx <- length(object[["L.f"]]) + 1
    L.f.flag[L.f.idx] <- TRUE
    thr.flag <- object[["E.v"]] <=
        (object[["E.v"]][L.v.idx] + object[["s.tr"]][L.v.idx])
    lm.L.v <- object[["lm.L.v"]]
    lm.L.f <- object[["lm.L.f"]]
    res <- list(summ.full = summary(object[["lm.full"]], ...),
                summ.L.v = if (!is.null(lm.L.v)) summary(lm.L.v, ...),
                summ.L.f = if (!is.null(lm.L.f)) summary(lm.L.f, ...),
                error.df = data.frame(n.inputs = 0:d,
                    E.tr = object[["E.tr"]], s.tr = object[["s.tr"]],
                    E.v = object[["E.v"]], L.f.flag = L.f.flag,
                    L.v.flag = L.v.flag, thr.flag = thr.flag))
    class(res) <- "summary.sisal"
    res
}

print.summary.sisal <- function(x, ...) {
    if(!inherits(x, "summary.sisal")) {
        stop('use only with "summary.sisal" objects')
    }
    underlined <- function(txt, symbol="-") {
        cat(txt, "\n", sep="")
        cat(paste0(rep(symbol, nchar(txt, type = "width") %/%
                       nchar(symbol, type = "width")), collapse=""),
            "\n", sep="")
    }

    underlined(gettext("summary() for linear models using all samples",
                       domain="R-sisal"), "=")

    underlined(gettext("All inputs", domain="R-sisal"))
    print(x[["summ.full"]], ...)

    underlined(gettext("L.v inputs", domain="R-sisal"))
    print(x[["summ.L.v"]], ...)

    underlined(gettext("L.f inputs", domain="R-sisal"))
    print(x[["summ.L.f"]], ...)

    error.df <- x[["error.df"]]
    multi.flag <- rowSums(cbind(error.df[["L.f.flag"]],
                                error.df[["L.v.flag"]],
                                error.df[["thr.flag"]]))
    error.df[["L.f.flag"]] <- NULL
    error.df[["L.v.flag"]] <- NULL
    error.df[["thr.flag"]] <- NULL
    error.df[["multi.flag"]] <-
        vapply(multi.flag, function(x) paste0(rep("*", x), collapse=""), "")
    names.e <- names(error.df)
    names.e[names.e == "n.inputs"] <- gettext("#Inputs", domain="R-sisal")
    names.e[names.e == "E.tr"] <- gettext("Tr. error", domain="R-sisal")
    names.e[names.e == "s.tr"] <- gettext("Tr. SD", domain="R-sisal")
    names.e[names.e == "E.v"] <- gettext("Val. error", domain="R-sisal")
    names.e[names.e == "multi.flag"] <- ""
    names(error.df) <- names.e

    underlined(gettext("Errors, measured with cross-validation",
                       domain="R-sisal"), "=")
    fArgs <- c(list(error.df),
               lapply(list(...), function(x) call("quote", x)))
    if (!("row.names" %in% names(fArgs))) {
        fArgs[["row.names"]] <- FALSE
    }
    do.call(print, fArgs)
    splitcat(gettext('Code: one "*" for each criterion:',
                     "* error is inside threshold",
                     "* smallest validation error (L.v)",
                     "* least complex model with error inside threshold (L.f)",
                     domain="R-sisal"))

    invisible(x)
}
