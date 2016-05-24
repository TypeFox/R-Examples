## Define check...() functions in a way they are compatible with same functions
## in the 'RUnit' package (these functions are directly inspired from those
## in RUnit). Make version that are more compatible with Komodo/SciViews-K Unit)

checkEquals <- function (target, current, msg = "",
tolerance = .Machine$double.eps^0.5, checkNames = TRUE, ...)
{
    val <- FALSE
    timing <- as.numeric(system.time({
        ret <- try({
            ## Run the test
            if (isTRUE(checkNames)) {
            	cn <- ""	# Since this is the default value
            } else {
            	cn <- ", checkNames = FALSE"
            	names(target) <- NULL
            	names(current) <- NULL
            }
            if (!is.numeric(tolerance))
                stop("tolerance has to be a numeric value")
            if (length(tolerance) != 1)
            	stop("tolerance has to be a scalar")
            res <- all.equal(target, current, tolerance = tolerance, ...)
            val <- isTRUE(res)
        }, silent = TRUE)
    }, gcFirst = FALSE)[3])
    ## Log this test
    test <- .logTest(timing)
    ## Decide if recording more info or not
    minTiming <- getOption("svUnit.minTiming")
    if (is.null(minTiming)) minTiming <- 0.1
    if (!isTRUE(getOption("svUnit.recordAll"))  && isTRUE(timing < minTiming)
        && val) return(invisible(TRUE))
    ## Check for error
    if (inherits(ret, "try-error")) {
        val <- NA
        .logTestData(test, msg = msg, call =
            deparse(sys.call()[1:3], nlines = 1), timing = timing, val = -1,
            res = as.character(ret))
    } else {
        .logTestData(test, msg = msg, call =
            deparse(sys.call()[1:3], nlines = 1), timing = timing, val = val,
            res = if (val) "" else paste(c(res, .formatResult(current)),
            collapse = "\n"))
    }
    return(invisible(val))
}

checkEqualsNumeric <- function (target, current, msg = "",
tolerance = .Machine$double.eps^0.5, ...)
{
    val <- FALSE
    timing <- as.numeric(system.time({
        ret <- try({
            ## Run the test
            if (!is.numeric(tolerance))
                stop("tolerance has to be a numeric value")
            if (length(tolerance) != 1)
                stop("tolerance has to be a scalar")
            res <- all.equal.numeric(as.vector(target), as.vector(current),
                tolerance = tolerance, ...)
            val <- isTRUE(res)
        }, silent = TRUE)
    }, gcFirst = FALSE)[3])
    ## Log this test
    test <- .logTest(timing)
    ## Decide if recording more info or not
    minTiming <- getOption("svUnit.minTiming")
    if (is.null(minTiming)) minTiming <- 0.1
    if (!isTRUE(getOption("svUnit.recordAll"))  && isTRUE(timing < minTiming)
        && val) return(invisible(TRUE))
    ## Check for error
    if (inherits(ret, "try-error")) {
        val <- NA
        .logTestData(test, msg = msg, call =
            deparse(sys.call()[1:3], nlines = 1), timing = timing, val = -1,
            res = as.character(ret))
    } else {
        .logTestData(test, msg = msg, call =
            deparse(sys.call()[1:3], nlines = 1), timing = timing, val = val,
            res = if (val) "" else paste(c(res, .formatResult(current)),
            collapse = "\n"))
    }
    return(invisible(val))
}

checkIdentical <- function (target, current, msg = "")
{
    val <- FALSE
    timing <- as.numeric(system.time({
        ret <- try({
            ## Run the test
            val <- identical(target, current)
        }, silent = TRUE)
    }, gcFirst = FALSE)[3])
    ## Log this test
    test <- .logTest(timing)
    ## Decide if recording more info or not
    minTiming <- getOption("svUnit.minTiming")
    if (is.null(minTiming)) minTiming <- 0.1
    if (!isTRUE(getOption("svUnit.recordAll"))  && isTRUE(timing < minTiming)
        && val) return(invisible(TRUE))
    ## Check for error
    if (inherits(ret, "try-error")) {
        val <- NA
        .logTestData(test, msg = msg, call =
            deparse(sys.call()[1:3], nlines = 1), timing = timing, val = -1,
            res = as.character(ret))
    } else {
        .logTestData(test, msg = msg, call =
            deparse(sys.call()[1:3], nlines = 1), timing = timing, val = val,
            res = .formatResult(current))
    }
    return(invisible(val))
}

checkTrue <- function (expr, msg = "")
{
    val <- FALSE
    timing <- as.numeric(system.time({
        ret <- try({
            ## Run the test
            val <- isTRUE(all(expr == TRUE))
        }, silent = TRUE)
    }, gcFirst = FALSE)[3])
    ## Log this test
    test <- .logTest(timing)
    ## Decide if recording more info or not
    minTiming <- getOption("svUnit.minTiming")
    if (is.null(minTiming)) minTiming <- 0.1
    if (!isTRUE(getOption("svUnit.recordAll"))  && isTRUE(timing < minTiming)
        && val) return(invisible(TRUE))
    ## Get call, without msg
    call <- sys.call()
    call <- deparse(call[names(call) != "msg"])
    ## Check for error
    if (inherits(ret, "try-error")) {
        val <- NA
        .logTestData(test, msg = msg, call =
            deparse(sys.call()[1:2], nlines = 1), timing = timing, val = -1,
            res = as.character(ret))
    } else {
        .logTestData(test, msg = msg, call =
        deparse(sys.call()[1:2], nlines = 1), timing = timing, val = val,
        res = .formatResult(expr))
    }
    return(invisible(val))
}

checkException <- function (expr, msg = "",
silent = getOption("svUnit.silentException"))
{
    val <- FALSE
    timing <- as.numeric(system.time({
        ret <- try({
            ## Run the test
            silent <- (is.null(silent) || isTRUE(silent))
            val <- inherits(res <- try(expr, silent = silent), "try-error")
        }, silent = TRUE)
    }, gcFirst = FALSE)[3])
    ## Log this test
    test <- .logTest(timing)
    ## Decide if recording more info or not
    minTiming <- getOption("svUnit.minTiming")
    if (is.null(minTiming)) minTiming <- 0.1
    if (!isTRUE(getOption("svUnit.recordAll"))  && isTRUE(timing < minTiming)
        && val) return(invisible(TRUE))
    ## Check for error
    if (inherits(ret, "try-error")) {
        val <- NA
        .logTestData(test, msg = msg, call =
            deparse(sys.call()[1:2], nlines = 1), timing = timing, val = -1,
            res = as.character(ret))
    } else {
        .logTestData(test, msg = msg, call =
            deparse(sys.call()[1:2], nlines = 1), timing = timing, val = val,
            res = if (val) paste(res, collapse = "\n") else
            "No exception generated!\n")
    }
    return(invisible(val))
}

DEACTIVATED <- function (msg = "")
    stop(msg)
