###################################################################
## cross-validation (bootstrap, k-fold cv etc.) of empirical risk
## for boosting algorithms for gamLSS models

make.grid <- function(max, length.out = 10, min = NULL, log = TRUE,
                      dense_mu_grid = TRUE) {

    if (is.null(min))
        min <- rep(1, length(max))
    if (length(min) == 1)
        min <- rep(min, length(max))
    if (length(length.out) == 1)
        length.out <- rep(length.out, length(max))
    if (length(length.out) != length(max))
        stop(sQuote("length.out"),
             " must be either scalar or a vector of the same length as ",
             sQuote("max"))
    if (length(min) != length(max))
        stop(sQuote("min"),
             " must be either scalar or a vector of the same length as ",
             sQuote("max"))

    if (log == TRUE) {
        min <- log(min)
        max <- log(max)
    }

    if (any(sapply(1:length(max), function(i) min[i] >= max[i])))
        stop("All min values must be smaller than the respectiv max value.")

    ## single paramter family
    if (length(max) == 1) {
        RET <- seq(from = min, to = max, length.out = length.out)
        if (log == TRUE)
            RET <- exp(RET)
        ## round to get integer values
        RET <- round(RET)
        if (any(duplicated(RET))) {
            warning("Duplicates produced; Only unique values are returned")
            RET <- unique(RET)
        }
        return(RET)
    }

    ## ELSE: multiple parameter family

    if (is.null(names(max)))
        stop(sQuote("max"), " must be a named vector")

    RET <- lapply(1:length(max), function(i)
                  seq(from = min[i], to = max[i],
                      length.out = length.out[i]))
    if (log == TRUE)
        RET <- lapply(RET, exp)
    ## round to get integer values
    RET <- lapply(RET, round)

    if (any(sapply(RET, function(x) any(duplicated(x))))) {
        warning("Duplicates produced; Only unique values are returned")
        RET <- lapply(RET, unique)
    }
    ## make grid
    RET <- expand.grid(RET)
    if (!is.null(names(max)))
        colnames(RET) <- names(max)

    ## no way to produce dense grid
    if (dense_mu_grid && max(RET[,1]) <= min(RET[, -1]))
        dense_mu_grid <- FALSE

    ## if dense grid for mu is requested:
    if (dense_mu_grid) {
        ## find all grid points where at least one mu grid point is greater than
        ## all other parameters
        tmp <- RET[RET[, 1] > apply(RET[, -1, drop = FALSE], 1, max), ]
        tmp <- unique(tmp[, -1, drop = FALSE])
        ## for these values produce a dense grid
        stop <- max(RET[, 1])
        for (i in 1:nrow(tmp)) {
            start <- max(tmp[i, ])
            res <- suppressWarnings(cbind(seq(start, to = stop, by = 1), tmp[i, ]))
            colnames(res) <- colnames(RET)
            RET <- rbind(RET, res)
        }
        RET <- unique(RET)
        ## SORT THE GRID AND THEN CONTINUE IN cvrisk BY checking if
        ##   grid[, 1] >= grid[, -1] && ???
        RET <- RET[do.call(order, RET[, rev(colnames(RET))]), ]
        # RET <- RET[order(RET[,2], RET[,1]), ]
        rownames(RET) <- NULL
    }
    attr(RET, "dense_mu_grid") <- dense_mu_grid
    return(RET)
}

###
# cvrisk, adapted version from mboost (2.2-2)
cvrisk.mboostLSS <- function(object, folds = cv(model.weights(object)),
                             grid = make.grid(mstop(object)),
                             papply = mclapply, trace = TRUE,
                             fun = NULL, ...) {

    weights <- model.weights(object)
    if (any(weights == 0))
        warning("Zero weights in ", sQuote("object"))
    if (is.null(folds)) {
        folds <- rmultinom(25, length(weights), weights/sum(weights))
    } else {
        stopifnot(is.matrix(folds) && nrow(folds) == length(weights))
    }
    if (length(object) != ncol(grid))
        stop(sQuote("grid"),
             " must be a matrix with the same number of columns",
             " as parameters in", sQuote("object"))

    if (!is.null(fun))
        stopifnot(is.function(fun))

    ### WHAT ABOUT:
    ## fam_name <- object$family@name
    call <- deparse(attr(object, "call"))
    oobrisk <- matrix(0, nrow = ncol(folds), ncol = ncol(grid))
    if (!is.null(attr(grid, "dense_mu_grid"))) {
        dense_mu_grid <- attr(grid, "dense_mu_grid")
    } else {
        dense_mu_grid <- FALSE
    }
    if (trace)
        cat("Starting cross-validation...\n",
            "[fold]\t[current mstop]\n", sep = "")
    if (is.null(fun)) {
        dummyfct <- function(i, weights, oobweights) {
            ## make model with new weights and minimal mstop
            mod <- update(object, weights = weights, oobweights = oobweights,
                          risk = "oobag", trace = FALSE,
                          mstop = apply(grid, 2, min))

            ## now we need to increase mstop (stupid or clever)
            risks <- vector("numeric", nrow(grid))

            ## fitting loop
            j <- 1
            while (j <= nrow(grid)) {
                j_start <- j
                ## check for dense grid
                if (dense_mu_grid &&
                    (grid[j, 1] >= max(grid[j, -1])) &&
                    j < nrow(grid) && ## needed for the next line (i.e., grid[j+1,])
                    (grid[j, 1] == grid[j + 1, 1] - 1)) {
                    ## now check how long the dense grid continues
                    # usual length:
                    j_tmp <- j + max(grid[, 1]) - grid[j, 1]
                    # now check if this is correct:
                    if (all(grid[j:j_tmp, -1] == grid[j, -1, drop = TRUE])) {
                        ## increase j
                        j <- j_tmp
                    }
                    ## else continue step by step...
                }
                mod[grid[j, ]]
                rsk <- risk(mod, merge = TRUE)
                risks[j_start:j] <- rsk[(length(rsk) - j + j_start):length(rsk)]
                if (trace) {
                    txt <- paste0(" [", i, "]\t",
                                  paste0("[", paste(mstop(mod), collapse = ","),
                                         "]"), "\n")
                    cat(txt)
                }
                j <- j + 1
            } ## end while
            return(risks)
        }
    }
    else {
        stop("currently not implemented")
        dummyfct <- function(i, weights, oobweights) {
            ## make model with new weights and minimal mstop
            mod <- update(object, weights = weights, oobweights = oobweights,
                          risk = "oobag", trace = FALSE,
                          mstop = apply(grid, 2, min))

            res <- vector("list", nrow(grid))
            for (i in 1:nrow(grid)) {
                mod[grid[i, ]]
                class(mod) <- class(object)
                res[i] <- fun(mod)
            }
            return(res)
        }
    }

    OOBweights <- matrix(rep(weights, ncol(folds)), ncol = ncol(folds))
    OOBweights[folds > 0] <- 0
    oobrisk <- papply(1:ncol(folds),
        function(i) dummyfct(i, weights = folds[, i],
                             oobweights = OOBweights[, i]), ...)
    ## get errors if mclapply is used
    if (any(idx <- sapply(oobrisk, is.character)))
        stop(sapply(oobrisk[idx], function(x) x))
    if (!is.null(fun))
        return(oobrisk)
    oobrisk <- t(as.data.frame(oobrisk))
    oobrisk <- oobrisk / colSums(OOBweights)
    colnames(oobrisk) <- apply(grid, 1,
                               function(x) paste(x, collapse = ","))
    rownames(oobrisk) <- 1:nrow(oobrisk)
    ## attr(oobrisk, "risk") <- fam_name
    ## sowas wie "Normal distribution: mu(id link)"
    attr(oobrisk, "call") <- call
    attr(oobrisk, "mstop") <- grid
    attr(oobrisk, "type") <- ifelse(!is.null(attr(folds, "type")),
        attr(folds, "type"), "user-defined")
    class(oobrisk) <- "cvriskLSS"
    oobrisk
}

print.cvriskLSS <- function(x, ...) {
    #cat("\n\t Cross-validated", attr(x, "risk"), "\n\t",
    cat("\n\t Cross-validated risk\n\t",
              attr(x, "call"), "\n\n")
    print(colMeans(x))
    cat("\n\t Optimal number of boosting iterations:", mstop(x), "\n")
    return(invisible(x))
}

plot.cvriskLSS <- function(x, type = c("heatmap", "lines"),
                           xlab = NULL, ylab = NULL,
                           ylim = range(x),
                           main = attr(x, "type"),
                            ...) {

    type <- match.arg(type)

    nms <- names(attr(x, "mstop"))
    if (type == "lines") {
        if (is.null(xlab))
            xlab <- paste0("Number of boosting iterations (",
                           paste0(nms, collapse = ","), ")")
        if (is.null(ylab))
            ylab <- "Out-of-bag risk"
    } else {
        if (is.null(xlab))
            xlab <- paste0("Number of boosting iterations (", nms[1], ")")
        if (is.null(ylab))
            ylab <- paste0("Number of boosting iterations (", nms[2], ")")
    }

    if (type == "lines") {
        cm <- colMeans(x)
        plot(1:ncol(x), cm, xlab = xlab, ylab = ylab,
             ylim = ylim, type = "n", lwd = 2,
             main = main, axes = FALSE, ...)
        out <- apply(x, 1, function(y) lines(1:ncol(x),y, col = "lightgrey"))
        rm(out)
        ms <- which.min(cm)
        lines(c(ms, ms), c(min(c(0, ylim[1] * ifelse(ylim[1] < 0, 2, 0.5))), cm[ms]),
              lty = 2)
        lines(1:ncol(x), cm, type = "l")
        axis(1, at = 1:ncol(x), labels = colnames(x))
        axis(2)
        box()
    } else {
        cm <- colMeans(x)
        grid <- attr(x, "mstop")
        if (!(ncol(grid) %in% c(2,3)))
            stop("currently only implemented for 2 and 3 dimensional grids")
        mstop <- mstop(x)
        #cm <- exp(cm)
        standardized_cm <- 1 - (cm - min(cm))/(max(cm) - min(cm))
        col <- grey(standardized_cm)
        # col <- heat.colors(length(unique(cm)))
        # col <- colorRampPalette(c("blue", "white", "red"))(length(cm))
        # plot(grid, col = col[order(cm)], pch = 15)
        # plot(grid, col = col, pch = 15)
        make_plot <- function(grid, col, idx, main, k) {
            plot(grid[idx, ], col = col[idx], xlab = xlab, ylab = ylab,
                 main = main, pch = 15, ...)
            if (is.null(k) || k == mstop[3]) {
                points(mstop[1], mstop[2], col = "red", pch = 22)
                lines(c(0, mstop[1]), c(mstop[2], mstop[2]), col = "red", lty = "dashed")
                lines(c(mstop[1], mstop[1]), c(0, mstop[2]), col = "red", lty = "dashed")
            }
        }
        if (ncol(grid) == 2) {
            make_plot(grid, col = col, idx = 1:nrow(grid), main = main, k = NULL)
        } else {
            for (k in unique(grid[, 3]))
                make_plot(grid[, -3], col = col, idx = (grid[, 3] == k),
                          main = paste0(main, "\n(", nms[3], "=", k, ")"),
                          k = k)
        }
    }
}

mstop.cvriskLSS <- function(object, parameter = NULL, ...) {
    res <- unlist(attr(object, "mstop")[which.min(colSums(object)),])
    if (!is.null(parameter)) {
        if(is.character(parameter))
            parameter <- extract_parameter(object, parameter)
        res <- res[parameter]
    }
    return(res)
}

if (FALSE) {
    library(gamboostLSS)

    ## check cvrisk
    set.seed(1907)
    x1 <- rnorm(1000)
    x2 <- rnorm(1000)
    x3 <- rnorm(1000)
    x4 <- rnorm(1000)
    x5 <- rnorm(1000)
    x6 <- rnorm(1000)
    mu    <- exp(1.5 +1 * x1 +0.5 * x2 -0.5 * x3 -1 * x4)
    sigma <- exp(-0.2 * x3)
    y <- numeric(1000)
    for( i in 1:1000)
        y[i] <- rnbinom(1, size = sigma[i], mu = mu[i])
    dat <- data.frame(x1, x2, x3, x4, x5, x6, y)
    model <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                         control = boost_control(mstop = 400),
                         center = TRUE)
    grid <- make.grid(c(mu = 1000, sigma = 1000), length.out = 5)
    plot(grid)
    abline(0,1)
    cvr <- cvrisk(model, folds = cv(model.weights(model), B = 5), grid = grid,
                  papply = lapply)

    ### check timings:

    ### Data generating process:
    set.seed(1907)
    x1 <- rnorm(1000)
    x2 <- rnorm(1000)
    x3 <- rnorm(1000)
    x4 <- rnorm(1000)
    x5 <- rnorm(1000)
    x6 <- rnorm(1000)
    mu    <- exp(1.5 +1 * x1 +0.5 * x2 -0.5 * x3 -1 * x4)
    sigma <- exp(-0.4 * x3 -0.2 * x4 +0.2 * x5 +0.4 * x6)
    y <- numeric(1000)
    for( i in 1:1000)
        y[i] <- rnbinom(1, size = sigma[i], mu = mu[i])
    dat <- data.frame(x1, x2, x3, x4, x5, x6, y)

    system.time({
        ## linear model with y ~ . for both components: 1 boosting iterations
        model <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                             control = boost_control(mstop = 1),
                             center = TRUE)
        for (i in 10:1000) {
            model[c(i, 10)]
        }
    })
    ## langsamer als:
    system.time({
        ## linear model with y ~ . for both components: 1 boosting iterations
        model <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                             control = boost_control(mstop = 1),
                             center = TRUE)
        model[c(1000, 10)]

    })


    ## what about the warnings in 3d?
}
