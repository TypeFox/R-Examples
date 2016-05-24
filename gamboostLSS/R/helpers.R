## helper functions

check <- function(what, what_char, names) {

    errormsg <- paste0(sQuote(what_char), " can be either a scalar, a (named) vector or a (named) list",
                      " of ", what_char, " values with same names as ",  sQuote("families"), "in ",
                      sQuote("boost_control"))

    if (is.list(what)) {
        if (is.null(names(what)) && length(what) == length(names))
            names(what) <- names
        if (!all(names(what) %in% names) ||
            length(unique(names(what))) != length(names))
            stop(errormsg)
        what <- what[names] ## sort in order of families
        what <- unlist(what)
    } else {
        if(length(what) != 1 && length(what) != length(names))
            stop(errormsg)
        if (length(what) == 1) {
            what <- rep(what, length(names))
            names(what) <- names
        } else {
            if (is.null(names(what)))
                names(what) <- names
            if (!all(names(what) %in% names))
                stop(errormsg)
            what <- what[names] ## sort in order of families
        }
    }

    return(what)
}

## extract data from a gamboostLSS model (used in plot_PI)
get_data <- function(x, which = NULL) {
    data <- attr(x, "data")
    if (length(data) != 0) {
        data <- data[, which, drop = FALSE]
    } else {
        data <- try(sapply(which, get, env = parent.frame(2)),
                    silent = TRUE)
        if (inherits(data, "try-error"))
                stop("No data set found.")
        data <- as.data.frame(data)
    }
    return(data)
}

## extract family name from a gamboostLSS model (used in plot_PI)
get_families_name <- function(x) {
    attr(attr(x, "families"), "name")
}

## make generic get_qfun function
get_qfun <- function(x, ...)
    UseMethod("get_qfun")

## extract family name from a gamboostLSS model (used in plot_PI)
get_qfun.mboostLSS <- function(x) {
    qfun <- attr(attr(x, "families"), "qfun")
    if (is.null(qfun))
        stop("Currently not implemented for this family")
    return(qfun)
}

## obtain pdf from gamlss.dist or global environment
## (needed in as.families)
get_qfun.character <- function(x) {
    qfun <- paste("gamlss.dist::q", x, sep = "")
    pdf <- try(eval(parse(text = qfun)), silent = TRUE)
    if (inherits(pdf, "try-error")) {
        ## try to find the function in global environment
        ## this is needed e.g. for truncated families
        qfun2 <- paste("q", x, sep = "")
        pdf <- try(eval(parse(text = qfun2)), silent = TRUE)
        if (inherits(pdf, "try-error"))
            stop(sQuote(qfun2), " and ", sQuote(qfun), " do not exist.")
    }
    return(pdf)
}

## obtain pdf from gamlss.dist or global environment
## (needed in as.families)
get_pdf <- function(x) {
    dfun <- paste("gamlss.dist::d", x, sep = "")
    pdf <- try(eval(parse(text = dfun)), silent = TRUE)
    if (inherits(pdf, "try-error")) {
        ## try to find the function in global environment
        ## this is needed e.g. for truncated families
        dfun2 <- paste("d", x, sep = "")
        pdf <- try(eval(parse(text = dfun2)), silent = TRUE)
        if (inherits(pdf, "try-error"))
            stop(sQuote(dfun2), " and ", sQuote(dfun), " do not exist.")
    }
    return(pdf)
}

## return mean or (first) modus of a vector depending on its class
mean_mod <- function(x) {
    if (is.numeric(x))
        return(mean(x, na.rm = TRUE))
    ## else compute and return modus
    if (is.character(x) || is.factor(x)) {
        ret <- names(which.max(table(x)))[1]
        if (is.factor(x))
            ret <- factor(ret, levels = levels(x))
        return(ret)
    }
    stop("not implemented for data type ", class(x))
}

## helper function copied from mboost_2.2-3
rescale_weights <- function(w) {
    if (max(abs(w - floor(w))) < sqrt(.Machine$double.eps))
        return(w)
    return(w / sum(w) * sum(w > 0))
}

## helper function in a modified version based on mboost_2.2-3
## print trace of boosting iterations
do_trace <- function(current, mstart, risk,
                     linebreak = options("width")$width / 2, mstop = 1000) {
    current <- current - mstart
    if (current != mstop) {
        if ((current - 1) %/% linebreak == (current - 1) / linebreak) {
            mchr <- formatC(current + mstart, format = "d",
                            width = nchar(mstop) + 1, big.mark = "'")
            cat(paste("[", mchr, "] ",sep = ""))
        } else {
            if ((current %/% linebreak != current / linebreak)) {
                cat(".")
            } else {
                cat(" -- risk:", risk[current + mstart], "\n")
            }
        }
    } else {
        cat("\nFinal risk:", risk[current + mstart], "\n")
    }
}

## helper function copied from mboost_2.2-3
### check measurement scale of response for some losses
check_y_family <- function(y, family)
    family@check_y(y)

################################################################################
# sapply function that differentiates between data.frames and (numeric) vectors
myApply <- function(X, FUN, ...) {
    ret <- lapply(X, FUN, ...)
    if (length(ret) == 1)
        ret <- ret[[1]]
    return(ret)
}


## helper function that stabilizes the negative gradient if requested by the user
stabilize_ngradient <- function(ngr, w = 1, stabilization) {
    ## set which to MAD if gamboostLSS_stab_ngrad = TRUE and which == "none"
    if (stabilization == "none" && getOption("gamboostLSS_stab_ngrad"))
        stabilization <- "MAD"
    ## stabilization using the mean absolute deviation (MAD)
    if (stabilization == "MAD") {
        div <- weighted.median(abs(ngr - weighted.median(ngr, w = w, na.rm = TRUE)),
                               w = w, na.rm = TRUE)
        div <- ifelse(div < 0.0001, 0.0001, div)
        ngr <- ngr / div
    }
    ngr
}


check_stabilization <- function(stabilization = c("none", "MAD")) {
    stabilization <- match.arg(stabilization)
    ## check if old stabilization interface is used and issue a warning
    if (getOption("gamboostLSS_stab_ngrad")) {
        warning("Usage of ", sQuote("options(gamboostLSS_stab_ngrad = TRUE)"),
                " is deprecated.\n", "Use argument ", sQuote("stabilization"),
                " in the fitting family. See ?Families for details.")
        if (stabilization == "none")
           warning(sQuote("stabilization"), " is set to ", dQuote("MAD"))
    }
    stabilization
}
