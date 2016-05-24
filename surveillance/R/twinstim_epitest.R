################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Monte Carlo Permutation Test for Space-Time Interaction in "twinstim"
###
### Copyright (C) 2015-2016 Sebastian Meyer
### $Revision: 1542 $
### $Date: 2016-02-01 15:10:18 +0100 (Mon, 01. Feb 2016) $
################################################################################

epitest <- function (model, data, tiles, method = "time", B = 199,
                     eps.s = NULL, eps.t = NULL, fixed = NULL,
                     verbose = TRUE, compress = FALSE, ...)
{
    ## check input
    stopifnot(inherits(model, "twinstim"), inherits(data, "epidataCS"),
              model$converged, isScalar(B), B >= 1)
    B <- as.integer(B)
    method <- match.arg(method, choices = c("LRT", "simulate", "time", "space"))
                                        # eval(formals(permute.epidataCS)$what)
    if (model$npars["q"] == 0L) {
        stop("no epidemic component in 'model'")
    }
    if (.epilink(model) == "log") {
        warning("boundary issues with the epidemic log-link",
                immediate. = TRUE)
    }
    if (isTRUE(fixed)) {
        fixed <- setdiff(grep("^e\\.", names(coef(model)), value = TRUE),
                         "e.(Intercept)")
    } else {
        stopifnot(is.null(fixed) || is.character(fixed))
    }
    t0 <- model$timeRange[1L]  # will not permute events before t0
    T <- model$timeRange[2L]
    
    ## auxiliary function to compute the LRT statistic
    lrt <- function (m0, m1) {
        l0 <- m0$loglik
        l1 <- m1$loglik
        c(l0 = l0, l1 = l1, D = 2 * (l1 - l0),
          converged = isTRUE(m1$converged) && isTRUE(m0$converged))
    }
    
    ## observed test statistic
    m0 <- update.twinstim(model, data = data,
                          epidemic = ~0, siaf = NULL, tiaf = NULL,
                          control.siaf = NULL, model = method == "simulate", cumCIF = FALSE,
                          cores = 1, verbose = FALSE,
                          optim.args = list(fixed = fixed, control = list(trace = 0)))
    if (!isTRUE(m0$converged)) {
        stop("endemic-only model did not converge")
    }
    LRT <- lrt(m0 = m0, m1 = model)
    STATISTIC_D <- structure(LRT["D"], l0 = LRT[["l0"]], l1 = LRT[["l1"]])
    STATISTIC_R0 <- c("simpleR0" = simpleR0(model, eps.s = eps.s, eps.t = eps.t))

    ## LRT p-value (CAVE: invalid for the default log-link models)
    DF <- length(coef(model)) - length(coef(m0)) # number of epidemic parameters
    PVAL_LRT <- pchisq(as.vector(STATISTIC_D), # drop attributes
                       df = DF, lower.tail = FALSE)
    
    ## result template
    res <- list(
        method = "Likelihood Ratio Test for Space-Time Interaction",
        data.name = paste0(deparse(substitute(data)),
            "\ntwinstim:  ", deparse(substitute(model))),
        statistic = STATISTIC_D,
        parameter = c("df" = DF),
        p.value = PVAL_LRT
    )
    class(res) <- c("epitest", "htest")

    if (method == "LRT") {
        ## we are done
        return(res)
    }
    
    ## otherwise: determine the null distribution via permutation or simulation
    res$method <- if (method == "simulate") {
        paste("Test for Space-Time Interaction (based on", B, "endemic simulations)")
    } else {
        "Monte Carlo Permutation Test for Space-Time Interaction"
    }
    if (model$npars["q"] > 1L) {
        warning("epidemic covariate effects might not be identifiable for null data",
                immediate. = TRUE)
    }

    ## define a function which generates data under the null
    generateNullData <- if (method == "simulate") {
        if (missing(tiles))
            stop("'tiles' is required for 'method = \"simulate\"'")
        rmarks <- .rmarks(data, t0 = t0, T = T)
        function() {
            events <- simEndemicEvents(m0, tiles = tiles)
            events@data <- cbind(events@data, rmarks(n = length(events)))
            as.epidataCS(events = events, stgrid = data$stgrid[,-1L],
                         W = data$W, qmatrix = data$qmatrix,
                         nCircle2Poly = attr(data$events$.influenceRegion, "nCircle2Poly"),
                         clipper = "polyclip", verbose = FALSE)
        }
    } else {
        function() permute.epidataCS(data, what = method, keep = time <= t0)
    }
    
    ## interpret 'verbose' level
    .verbose <- if (is.numeric(verbose)) {
        if (verbose >= 2) {
            ## create '.verbose' expression to print test statistics
            stats2string <- function (lrt, simpleR0)
                paste0(c(names(lrt)[1:3], "simpleR0"), " = ",
                       sprintf(paste0("%4.", c(0,0,1,2), "f"), c(lrt[1:3], simpleR0)),
                       collapse = " | ")
            cat("Endemic/Epidemic log-likelihoods, LRT statistic, and simple R0:\n",
                stats2string(LRT, STATISTIC_R0), "\n",
                "\nResults from B=", B, if (method == "simulate")
                    " endemic simulations" else paste0(" permutations of ", method),
                ## will actually not be printed if parallelized using clusters ...
                ":\n", sep = "")
            substitute({
                cat(STATS2STRING)
                if (!lrt["converged"]) {
                    msg <- c(m0 = m0$converged, m1 = m1$converged)
                    msg <- msg[msg != "TRUE"]
                    cat(" | WARNING (", paste0(names(msg), collapse = " and "),
                        "): ", paste0(unique(msg), collapse = " and "), sep = "")
                }
                cat("\n")
            }, list(STATS2STRING = body(stats2string)))
        } else {
            verbose <- verbose == 1
        }
    } else verbose

    siafInt <- NULL
    if (method != "simulate") {
        ## if siafpars are fixed, determine siafInt for use in all permutations
        siafpars <- coeflist(model)$siaf
        if (length(siafpars) > 0L && all(names(siafpars) %in% fixed) &&
            is.null(siafInt <- environment(model)$siafInt)) {
            if (!identical(FALSE, verbose))
                cat("pre-evaluating 'siaf' integrals with fixed parameters ...\n")
            setup <- update.twinstim(model, data = data, optim.args = NULL, verbose = FALSE)
            assign("siafpars", siafpars, envir = environment(setup))
            siafInt <- with(environment(setup), do.call("..siafInt", .siafInt.args))
        }
    }
    
    ## define the function to be replicated B times:
    ## permute/simulate data, update epidemic model, compute endemic-only model,
    ## and compute test statistics
    permfits1 <- function (...) {
        ## depends on 'data', 'model', 'lrt', 'eps.s', 'eps.t', and 'fixed'
        .permdata <- generateNullData()
        .siafInt <- if (!is.null(siafInt)) {
            siafInt[match(row.names(.permdata$events), row.names(data$events))]
        } # else NULL
        ## sink(paste0("/tmp/trace_", Sys.getpid()), append = TRUE)
        m1 <- update.twinstim(model, data = .permdata,
                              control.siaf = list(siafInt = .siafInt),
                              model = FALSE, cumCIF = FALSE,
                              cores = 1, verbose = FALSE,
                              optim.args = list(fixed = fixed,
                                  control = list(trace = is.numeric(verbose) && verbose >= 3)))
        ## sink()
        m0 <- update.twinstim(m1, epidemic = ~0, siaf = NULL, tiaf = NULL,
                              control.siaf = NULL,
                              optim.args = list(control = list(trace = 0)))
        lrt <- lrt(m0, m1)
        simpleR0 <- simpleR0(m1, eps.s = eps.s, eps.t = eps.t)
        if (isTRUE(compress)) { # save memory
            m0[c("fitted", "fittedComponents", "R0")] <-
                m1[c("fitted", "fittedComponents", "R0")] <- list(NULL)
        }
        list(m0 = m0, m1 = m1,
             stats = c(lrt[1:3], simpleR0 = simpleR0, lrt["converged"]))
    }

    ## rock'n'roll (the computationally intensive part)
    permfits <- plapply(X = integer(B), FUN = permfits1,
                        .verbose = .verbose, ...)

    ## if parallelized using forking with insufficient memory available,
    ## part of the replications in 'permfits' may be left unassigned (NULL)
    permIsNull <- vapply(X = permfits, FUN = is.null,
                         FUN.VALUE = logical(1L), USE.NAMES = FALSE)
    if (npermIsNull <- sum(permIsNull)) {
        warning(npermIsNull, "/", B,
                " replications did not return (insufficient memory?)")
        permfits <- permfits[!permIsNull]
    }
    
    ## extract the statistics
    permstats <- as.data.frame(t(vapply(
        X = permfits, FUN = "[[", "stats",
        FUN.VALUE = numeric(5L), USE.NAMES = TRUE
    )))
    permstats$converged <- as.logical(permstats$converged)
    
    ## compute permutation-based p-value
    PVAL_D <- mean(c(STATISTIC_D, permstats[permstats$converged, "D"]) >= STATISTIC_D)
    PVAL_R0 <- mean(c(STATISTIC_R0, permstats[permstats$converged, "simpleR0"]) >= STATISTIC_R0)
    
    ## set results
    res$statistic <- structure(STATISTIC_R0, "D" = unname(STATISTIC_D))
    res$parameter <- c("B" = sum(permstats$converged))
    res$p.value <- structure(PVAL_R0, "D-based" = PVAL_D, "LRT" = PVAL_LRT)
    res$permfits <- permfits
    res$permstats <- permstats
    res
}

coef.epitest <- function (object, which = c("m1", "m0"), ...)
{
    which <- match.arg(which)
    permcoefs <- vapply(X = object$permfits,
                        FUN = function (x) coef(x[[which]]),
                        FUN.VALUE = coef(object$permfits[[1L]][[which]]),
                        USE.NAMES = TRUE)
    t(permcoefs)
}

plot.epitest <- function (x, teststat = c("simpleR0", "D"), ...)
{
    teststat <- match.arg(teststat)
    defaultArgs <- switch(teststat,
        "simpleR0" = list(
            permstats = x$permstats$simpleR0,
            xmarks = setNames(x$statistic, "observed"),
            xlab = expression("Simple " * R[0])
        ),
        "D" = list(
            permstats = x$permstats$D,
            xmarks = setNames(attr(x$statistic, "D"), "observed"),
            xlab = expression(D == 2 %.% log(L[full]/L[endemic]))
        )
    )
    args <- modifyList(defaultArgs, list(...))
    if (is.null(args[["permstats"]]))
        stop("nothing to plot (no 'permstats' available)")
    do.call("permtestplot", args)
}

## auxiliary function also used by plot.knox(), permutationTest(), ...
permtestplot <- function (permstats, xmarks = NULL, xlab = "test statistic", ...)
{
    defaultArgs <- list(
        data = permstats, xlab = xlab, col = "lavender",
        main = "Monte Carlo permutation test for space-time interaction",
        xlim = extendrange(c(permstats, xmarks))
    )
    do.call("truehist", modifyList(defaultArgs, list(...), keep.null = TRUE))
    if (!is.null(xmarks)) {
        abline(v = xmarks, lwd = 2)
        axis(3, at = xmarks, labels = names(xmarks), # if NULL the value is used
             tick = FALSE, line = -1, font = 2)
    }
    invisible(NULL)
}
