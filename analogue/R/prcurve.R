## Fit a Principal Curve to matrix X
## Wrapper to principal.curve() in package princurve
## We use the original code plus our wrappers as pcurve()
## in package pcurve is too complex for our needs

## prcurve (named after prcomp): fits a principal curve to matrix X
prcurve <- function(X,
                    method = c("ca","pca","random","user"),
                    start = NULL,
                    smoother = smoothSpline,
                    complexity,
                    vary = FALSE,
                    maxComp,
                    finalCV = FALSE,
                    axis = 1,
                    rank = FALSE,
                    stretch = 2,
                    maxit = 10,
                    trace = FALSE,
                    thresh = 0.001,
                    plotit = FALSE,
                    ## fitFUN = c("princurve","pcurve"),
                    ## latent = FALSE,
                    ...) {
    ## X should be a matrix, attempt to coerce
    if(!isTRUE(inherits(X, "matrix")))
        X <- data.matrix(X)
    ## set/select default method for starting configuration
    if(missing(method))
        method <- "ca"
    else
        method <- match.arg(method)
    ## ## set/select default fitting function
    ## if(missing(fitFUN))
    ##     fitFUN <- "princurve"
    ## else
    ##     fitFUN <- match.arg(fitFUN)
    ## if(latent && fitFUN == "princurve")
    ##     warning("Scaling PC to a latent variable not availble with fitFUN = \"princurve\".")
    ## data stats
    n <- NROW(X) ## number of observations
    m <- NCOL(X) ## number of variables

    ## fit a PCA and store in result
    ord <- rda(X)

    ## starting configuration
    config <- startConfig <- initCurve(X, method = method,
                                       rank = rank,
                                       axis = axis,
                                       start = start)
    ## Need to sort out auto DF choices after pcurve::pcurve
    ## Vary degrees of freedom per variable?
    if(missing(complexity)) {
        complexity <- numeric(length = m)
        if (trace){ ## set up progress bar
            writeLines("\n   Determining initial DFs for each variable...")
            pb <- txtProgressBar(max = m, style = 3)
        }
        for(j in seq_along(complexity)) {
            if(trace) { ## update progress
                setTxtProgressBar(pb, j)
            }
            ## fit the mode & grab DF
            complexity[j] <-
                smoother(config$lambda, X[, j], choose = TRUE, ...)$complexity
        }
        if (trace) { ## finalise the progress bar
            close(pb)
            writeLines("\n")
        }

        if(!vary) { ## median complexity for all vars
            complexity <- rep(median(complexity), m)
        }
    } else {
        if((len <- length(complexity)) == 1) {
            complexity <- rep(complexity, m)
        } else if(len != m) {
            stop("Ambiguous 'complexity'; should be length 1 or NCOL(X)")
        }
    }
    if(missing(maxComp))
        maxComp <- 5 * log10(n)
    ## fix-up/reset complexity > maxComp to maxComp
    complexity[complexity > maxComp] <- maxComp
    ##
    iter <- 0L
    if(trace) {
        ##writeLines(strwrap(tmp <- paste(rep("-", options("width")[[1]]),
        ##                                collapse = "")))
        writeLines("Fitting Principal Curve:\n")
        writeLines(sprintf("Initial curve: d.sq: %.3f", config$dist))
    }

    ## vary == FALSE needs to set some things for smoothers like GAM
    ## which will select smoothness even if complexity stated
    smooths <- c("smoothGAM")
    if(!vary && (deparse(substitute(smoother)) %in% smooths)) {
        CHOOSE <- TRUE
    } else {
        CHOOSE <- FALSE
    }

    ##dist.raw <- sum(diag(var(X))) * (NROW(X) - 1)
    dist.old <- sum(diag(var(X)))
    s <- matrix(NA, nrow = n, ncol = m)
    converged <- (abs((dist.old - config$dist)/dist.old) <=
                  thresh)
    ## Start iterations ----------------------------------------------
    ### - store fitted smoothers in list
    smooths <- vector(mode = "list", length = m)
    while (!converged && iter < maxit) {
        iter <- iter + 1L
        for(j in seq_len(m)) {
            smooths[[j]] <- smoother(config$lambda, X[, j],
                                     complexity = complexity[j],
                                     choose = CHOOSE, ...)
            s[, j] <- fitted(smooths[[j]])
        }
        ##
        dist.old <- config$dist
        ## if(fitFUN == "princurve") {
            config <- get.lam(X, s = s, stretch = stretch)
        ## } else {
        ##     uni.lam <- sort(unique(config$lambda))
        ##     config <- pcget.lam(X, s = s, latent = latent, stretch = stretch,
        ##                         uni.lam = uni.lam)
        ## }
        class(config) <- "prcurve"
        ## Converged?
        converged <- (abs((dist.old - config$dist)/dist.old) <=
                      thresh)
        if(plotit) {
            ## plot the iteration -- need to add some components
            ## because of changes to plot method
            dev.hold()
            config$data <- X
            config$ordination <- ord
            plot(config, sub = paste("Iteration:", iter))
            dev.flush()
        }
        if (trace)
            writeLines(sprintf(paste("Iteration %",
                                     max(3, nchar(maxit)),
                                     "i: d.sq: %.3f", sep = ""),
                               iter, config$dist))
    }
    ## End iterations ------------------------------------------------
    ## if we want a final CV spline fit?
    if(finalCV) {
        iter <- iter + 1L
        for(j in seq_len(n)) {
          smooths[[j]] <- smoother(config$lambda, X[, j],
                                   cv = TRUE, choose = TRUE, ...)
          if(smooths[[j]]$complexity > maxComp) {
            smooths[[j]] <- smoother(config$lambda, X[, j], cv = FALSE,
                                     choose = FALSE,
                                     complexity = maxComp,
                                     ...)
          }
          s[, j] <- fitted(smooths[[j]])
            ## sFit <- smoother(config$lambda, X[, j],
            ##                  cv = TRUE, choose = TRUE, ...)
            ## s[, j] <- if(sFit$complexity > maxComp) {
            ##     ## too complex, turn of CV and refit with max df allowed
            ##     fitted(smoother(config$lambda, X[, j], cv = FALSE,
            ##                     choose = FALSE,
            ##                     complexity = maxComp,
            ##                     ...))
            ## } else {
            ##     fitted(sFit)
            ## }
        }
        config <- get.lam(X, s = s, stretch = stretch)
        class(config) <- "prcurve"
        if(plotit) {
            ## plot the iteration -- need to add some components
            ## because of changes to plot method
            dev.hold()
            config$data <- X
            config$ordination <- ord
            plot(config)
            dev.flush()
        }
        if (trace)
            writeLines(sprintf(paste("Iteration %", max(3, nchar(maxit)),
                                     "s: d.sq: %.4f", sep = ""),
                               "CV", config$dist))
    }
    if(trace){
        cat("\n")
        if(converged) {
            writeLines(strwrap(paste("PC Converged in", iter, "iterations.")))
        } else {
            writeLines(strwrap(paste("PC did not converge after", iter,
                                     "iterations.")))
        }
        cat("\n")
    }
    ## prepare objects for return
    names(config$tag) <- names(config$lambda) <-
        rownames(config$s) <- rownames(X)
    colnames(config$s) <- names(complexity) <- colnames(X)
    config$converged <- converged
    config$iter <- iter
    config$totalDist <- startConfig$dist
    config$complexity <- complexity
    ## config$fitFUN <- fitFUN
    config$smooths <- smooths
    names(config$smooths) <- colnames(X)
    config$call <- match.call()
    config$ordination <- ord
    config$data <- X
    config$stretch <- stretch
    class(config) <- c("prcurve")
    config
}

`print.prcurve` <- function(x, digits = max(3, getOption("digits") - 3),
                            ...) {
    cat("\n")
    writeLines(strwrap("Principal Curve Fitting", prefix = "\t"))
    cat("\n")
    writeLines(strwrap(pasteCall(x$call)))
    cat("\n")
    writeLines(strwrap(paste("Algorithm",
                             ifelse(x$converged, "converged", "failed to converge"),
                             "after",
                             x$iter,
                             ifelse(isTRUE(all.equal(x$iter, 1)),
                                    "iteration", "iterations"),
                             sep = " ")))
    cat("\n")
    tab <- cbind(c(x$totalDist, x$totalDist - x$dist, x$dist),
                 c(1.00, (x$totalDist - x$dist) / x$totalDist,
                   x$dist / x$totalDis))
    dimnames(tab) <- list(c("Total","Explained","Residual"),
                          c("SumSq","Proportion"))
    printCoefmat(tab, digits = digits, na.print = "")
    cat("\n")
    writeLines(strwrap(paste("Fitted curve uses",
                             round(edf <- sum(x$complexity), digits = digits),
                             ifelse(edf > 1, "degrees", "degree"),
                             "of freedom.", sep  = " ")))
    invisible(x)
}
