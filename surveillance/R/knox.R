################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Knox test for space-time interaction
###
### Copyright (C) 2015 Sebastian Meyer
### $Revision: 1347 $
### $Date: 2015-05-29 11:45:51 +0200 (Fre, 29. Mai 2015) $
################################################################################


knox <- function (dt, ds, eps.t, eps.s, simulate.p.value = TRUE, B = 999, ...)
{
    stopifnot(length(dt) == length(ds))

    ## logical vectors indicating which pairs are close in time and space
    closeInTime <- if (is.logical(dt)) {
        dt
    } else {
        stopifnot(is.numeric(dt), isScalar(eps.t))
        dt <= eps.t
    }
    closeInSpace <- if (is.logical(ds)) {
        ds
    } else {
        stopifnot(is.numeric(ds), isScalar(eps.s))
        ds <= eps.s
    }

    ## manually build the contingency table (table() with factor() is too slow)
    .lab <- c("close", "not close")
    knoxtab <- array(
        tabulate(4L - closeInTime - 2L*closeInSpace, nbins = 4L),
        dim = c(2L, 2L),
        dimnames = list(
            dt = if (is.logical(dt)) .lab else paste(c("<=", " >"), eps.t),
            ds = if (is.logical(ds)) .lab else paste(c("<=", " >"), eps.s)
        ))
    class(knoxtab) <- "table"

    ## expected number of close pairs in the absence of spatio-temporal interaction
    npairs <- sum(knoxtab)
    expected <- sum(knoxtab[1L,]) / npairs * sum(knoxtab[,1L])
    ##<- this order of terms avoids integer overflow
    
    ## test statistic is the number of spatio-temporally close pairs
    METHOD <- "Knox test"
    STATISTIC <- knoxtab[1L]

    ## determine statistical significance
    pval_Poisson <- ppois(STATISTIC, expected, lower.tail = FALSE)
    PVAL <- if (simulate.p.value) { # Monte Carlo permutation approach
        stopifnot(isScalar(B))
        B <- as.integer(B)
        METHOD <- paste(METHOD, "with simulated p-value")
        PARAMETER <- setNames(B, "B")
        permstats <- plapply(X = integer(B), FUN = function (...)
            sum(closeInSpace & closeInTime[sample.int(npairs)]), ...)
        structure(mean(c(STATISTIC, permstats, recursive = TRUE) >= STATISTIC),
                  Poisson = pval_Poisson)
    } else {
        METHOD <- paste(METHOD, "with Poisson approximation")
        PARAMETER <- setNames(expected, "lambda")
        pval_Poisson
    }

    ## return test results
    structure(
        list(method = METHOD,
             data.name = paste("dt =", deparse(substitute(dt)),
                               "and ds =", deparse(substitute(ds))),
             statistic = setNames(STATISTIC, "number of close pairs"),
             parameter = PARAMETER, p.value = PVAL, alternative = "greater",
             null.value = setNames(expected, "number"),
             permstats = if (simulate.p.value) {
                 unlist(permstats, recursive = FALSE, use.names = FALSE)
             },
             table = knoxtab),
        class = c("knox", "htest")
    )
}

print.knox <- function (x, ...)
{
    ## first print by the default method for class "htest"
    NextMethod("print")

    ## then also output the contingency table
    cat("contingency table:\n")
    print(x$table)
    cat("\n")
    invisible(x)
}

plot.knox <- function (x, ...)
{
    if (is.null(permstats <- x[["permstats"]])) {
        stop("this plot-method is for a permutation-based Knox test")
    }
    defaultArgs <- list(
        permstats = permstats,
        xmarks = setNames(c(x[["null.value"]], x[["statistic"]]),
            c("expected", "observed")),
        xlab = "number of close pairs"
    )
    do.call("permtestplot", modifyList(defaultArgs, list(...)))
}

xtable.knox <- function (x, caption = NULL, label = NULL,
                         align = paste0("r|rr", if (!is.null(sumlabel)) "|r"),
                         digits = 0, display = NULL, ...,
                         sumlabel = "$\\sum$")
{
    tab <- x$table
    if (!is.null(sumlabel)) {
        FUN <- setNames(list(sum), sumlabel)
        tab <- addmargins(tab, FUN = FUN, quiet = TRUE)
    }
    xtable(tab, caption = caption, label = label, align = align,
           digits = digits, display = display, ...)
}

toLatex.knox <- function (object, dnn = names(dimnames(object$table)),
                          hline.after = NULL, sanitize.text.function = NULL, ...)
{
    xtab <- xtable(object, ...)
    if (is.null(hline.after))
        hline.after <- unique(c(-1,0,2,nrow(xtab)))
    if (is.null(sanitize.text.function))
        sanitize.text.function <- function (x)
            gsub("<=", "$\\le$", gsub(">", "$>$", x, fixed = TRUE), fixed = TRUE)
    res <- toLatex.xtable(xtab, hline.after = hline.after,
                          sanitize.text.function = sanitize.text.function, ...)
    if (is.null(dnn)) {
        res
    } else {
        stopifnot(length(dnn) == 2)
        headeridx <- grep("&", res, fixed = TRUE)[1L]
        res[headeridx] <- paste0(dnn[1L], res[headeridx])
        res <- append(res, paste0(" & \\multicolumn{2}{|c|}{", dnn[2L], "} & \\\\"),
                      after = headeridx - 1L)
        class(res) <- "Latex"
        res
    }
}
