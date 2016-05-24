# model.matrix.earth.R: functions for manipulating earth model matrices
#
# The main functions are:
#
# expand.arg(x, env, is.y.arg) expand factors in x and convert to double mat with col names
#      Called by earth.formula, earth.default, get.earth.x
#
# model.matrix(terms, data) standard R function to expand factors
#     Called by expand.arg earth.formula, get.earth.x, predict.earth(type="terms")
#
# model.matrix.earth(object, x, ...) x arg must not be expanded, returns bx
#     Called by predict.earth
#
# get.earth.x(object, data) returns x expanded for factors and all double
#      Called by model.matrix.earth
#
# get.bx(x, which.terms, dirs, cuts) x arg must be already expanded
#      Called by model.matrix.earth, pruning.pass
#
#-----------------------------------------------------------------------------

# If model.frame can't interpret the data passed to it it silently
# returns the fitted values.  This routine makes that not silent.
# Note: model.frame is a standard R library function (it's in stats/R/models.R).

check.nrows <- function(expected.nrows, actual.nrows, fitted.nrows, Callers.name)
{
    if(actual.nrows != expected.nrows) {
        if(actual.nrows == fitted.nrows)
            stop0("model.frame.default could not interpret the data passed to ",
                  Callers.name,
                  "\n        (actual.nrows=", actual.nrows,
                  " expected.nrows=", expected.nrows,
                  " fitted.nrows=", fitted.nrows, ")")
        else  # can probably never get here
            warning0(Callers.name, " returned a number ", actual.nrows,
                     " of rows that was different from the number ",
                     expected.nrows,
                     " of rows in the data")
    }
}
good.colname <- function(name)
{
    # The nchar check prevents super long names (60 is arb)
    # that are actually contents of vectors e.g. c(1,2,3,etc.)
    # The grep ensures that there are no more than three commas,
    # also to prevent using the contents of vectors.

    !is.null(name) && nchar(name) <= 60 && !grepany(",.*,.*,", name)
}
good.colnames <- function(x)
{
    colnames <- colnames(x)
    if(is.null(colnames))
        return(FALSE)
    for(i in seq_along(colnames))
        if(!good.colname(colnames[i]))
            return(FALSE)
    return(TRUE)
}
# Called from earth.fit just before doing the pruning pass
# Also called by model.matrix.earth (which returns bx)
# The x arg must be already expanded

get.bx <- function(x, which.terms, dirs, cuts)
{
    stopifnot(all(dirs[1,] == 0))   # intercept term dirs must all be 0
    check.which.terms(dirs, which.terms)
    stopifnot(NCOL(x) > 0)
    colnames <- rownames(dirs[which.terms,,drop=FALSE])
    bx <- matrix(0, nrow=nrow(x), ncol=length(which.terms),
                 dimnames=list(NULL, colnames))
    ibx <- 1
    for(iterm in which.terms) {
        temp1 <- 1
        for(ipred in seq_len(ncol(x))) {
            dir <- dirs[iterm, ipred]
            if(dir == 2)  # predictor enters linearly?
                temp1 <- temp1 * x[, ipred]
            else if(dir == -1 || dir == 1) {
                temp2 <- dir * (x[, ipred] - cuts[iterm, ipred])
                temp1 <- temp1 * temp2 * (temp2 > 0)
            } else if(dir != 0)
                stop0("illegal direction ", dir, " in 'dirs'")
        }
        bx[, ibx] <- temp1
        ibx <- ibx + 1
    }
    bx
}
# Called only by model.matrix.earth

get.earth.x <- function(    # returns x expanded for factors
    object  = stop("no 'object' argument"),
    data    = NULL,         # can be a dataframe, matrix, or vector
    env,                    # environment for evaluation
    trace   = 0,
    Callers.name)           # caller's name for trace messages
{
    # Given an x mat, return an x mat with column names equal to
    # expected.colnames and with the columns in their correct order
    # So the user can hand us an x mat without column names,
    # or with named columns but in the wrong order, or an xmat containing
    # only a needed subset of all the columns, etc.
    # This code is a mess and doesn't handle all cases, just the common ones.

    fix.x.columns <- function(x, expected.colnames)
    {
        colnames <- colnames(x)
        ncolnames <- length(colnames)
        nexpected <- length(expected.colnames)
        if(is.null(colnames)) {
            if(trace >= 1)
                cat0(Callers.name, ": x has no column names, ",
                     "adding column names: ", paste.collapse(expected.colnames),
                     "\n")
            colnames <- expected.colnames
         } else if(ncolnames < nexpected) {
            # CHANGED Oct 2008: allow user to specify less than the expected
            # nbr of columns -- which is ok if he specifies all predictors
            # actually used by the model.

            imatch <- pmatch(colnames, expected.colnames, nomatch=0)
            if(any(imatch == 0)) {
                # can't repair the error because there are colnames in x that aren't
                # in expected.names (tends to happen with expanded factor names)
                stop0(Callers.name, ": x has ", ncolnames,
                      " columns, expected ", length(expected.colnames),
                      " to match: ", paste.collapse(expected.colnames))
            }
            # Create a new x, putting the existing cols into their correct positions.
            # Cols that aren't in the original x will end up as all 999s in the
            # the recreated x; that doesn't matter if they are for predictors
            # that are unused in the earth model.
            # We use 999 instead of NA else the model matrix routines fail incorrectly
            # because they don't like NAs when doing predictions.
            if(trace >= 1)
                cat0(Callers.name, ": x has missing columns, ",
                     "creating a new x with all cols\n")
            imatch <- pmatch(expected.colnames, colnames, nomatch=0)
            x.original <- x
            x <- matrix(data=999, nrow=nrow(x), ncol=nexpected)
            for(i in seq_len(nexpected))
                if(imatch[i])
                    x[,i] <- x.original[,imatch[i]]
            colnames <- expected.colnames
        } else if(ncolnames > nexpected)
            NULL # TODO not sure what to do here (do nothing so old regression tests pass)
        else {
            imatch <- pmatch(colnames, expected.colnames, nomatch=0)
            if(all(imatch == 0)) {
                   if(trace >= 1)
                        cat0(Callers.name,
                             ": unexpected x column names, renaming columns\n",
                             "    Old names: ", paste.collapse(colnames), "\n",
                             "    New names: ", paste.collapse(expected.colnames),
                             "\n")

                    colnames <- expected.colnames
            } else {
                 # replace indices for non-found predictor names with their value
                 # i.e. assume columns with unknown names are in their right position
                 for(i in seq_along(imatch))
                     if(imatch[i] == 0)
                         imatch[i] = i

                 # if any columns are in the wrong order then fix their order
                 # (imatch will be 1,2,3,... if columns are in the right order)

                 if(!all(imatch == seq_along(imatch))) {
                   s <- paste0(Callers.name, ": x columns are in the wrong order%s\n",
                               "    Old columns: ", paste.collapse(colnames), "\n",
                               "    New columns: ", paste.collapse(expected.colnames),
                               "\n")
                   if(length(imatch) == ncol(x)) {
                       trace1(trace, s, ", correcting the column order")
                       x <- x[,imatch]
                       colnames <- colnames[imatch]
                   } else
                       warnf(s, "")
                }
            }
        }
        colnames(x) <- colnames
        x
    }
    # Return x with matrix dimensions.
    # If x is already a matrix this does nothing.
    # Else allow x to be a vector if its length is an integer multiple
    # of the number of columns in the original x

    my.as.matrix <- function(x)
    {
        if(is.null(ncol(x))) {
            nrows <- length(x) / length(object$namesx)
            if(floor(nrows) == nrows)
                dim(x) <- c(nrow=nrows, ncol=length(x) / nrows)
            else
                stop0(Callers.name, ":\n",
                      "       could not convert vector x to matrix because ",
                      "length(x) ", length(x), "\n",
                      "       is not a multiple of the number ",
                      if(length(ncol(object$namesx))) ncol(object$namesx) else 0,
                      " of predictors ",
                      "\n       Expected predictors: ",
                      if(is.null(colnames(object$namesx)))
                          "none?"
                      else
                          paste.collapse(colnames(object$namesx))
                      )
        }
        x
    }
    check.expanded.ncols <- function(x, object)
    {
        if(NCOL(x) != NCOL(object$dirs)) {
            stop0(Callers.name,
                     ": the number ", NCOL(x), " of columns of x\n",
                     "(after factor expansion) does not match the number ",
                     NCOL(object$dirs),
                     " of columns of the earth object",
                     "\n    expanded x:  ", paste.collapse(colnames(x)),
                     "\n    object$dirs: ", paste.collapse(colnames(object$dirs)),
                     "\nPossible remedy: check factors in the input data")
        }
    }
    #--- get.earth.x starts here ---

    trace <- get.update.arg(trace, "trace", object, env,
                            trace1=NULL, Callers.name, print.trace=FALSE)
    if(is.null(trace))
        trace <- 0
    this.call <- match.call()
    if(is.null(object$terms)) {
        # object was created with earth.default, no formula

        x <- get.update.arg(data, "x", object, env, trace, Callers.name)
        x <- my.as.matrix(x)
        x <- fix.x.columns(x, object$namesx)
        if(trace >= 1) {
            print_summary(x, sprintf("%s: x", Callers.name), trace=2)
            trace2(trace, "\n")
        }
        x <- expand.arg(x, env, trace, is.y.arg=FALSE)
    } else {
        # object was created with earth.formula

        Terms <- delete.response(object$terms)
        data <- get.update.arg(data, "data", object, env, trace, Callers.name)
        data <- my.as.matrix(data)
        data <- fix.x.columns(data, object$namesx)
        data <- as.data.frame(data)
        expected.nrows <- nrow(data)
        if(trace >= 1) {
            print_summary(data, sprintf("%s: x", Callers.name), trace=2)
            trace2(trace, "\n")
        }
        data <- model.frame(Terms, data, na.action=na.pass)
        if(trace >= 1) {
            print_summary(data, sprintf("%s: after call to model.frame: data", Callers.name),
                          trace=2)
            trace2(trace, "\n")
        }
        classes <- attr(Terms, "dataClasses")
        if(!is.null(classes)) {
            # Use "try" for leniency, to allow numeric to be used for factors etc.
            # There is special treatment for the following message because it seems to be benign:
            #   variable 'foo' was fitted with type "nmatrix.1" but type "numeric" was supplied
            try <- try(.checkMFClasses(classes, data), silent=TRUE)
            if(is.try.err(try) && !grepl("\"nmatrix.1\" .* \"numeric\"", try[1])) {
                cat(try)
                cat("Continuing anyway, first few rows of x are\n")
                print(head(data))
            }
        }
        x <- model.matrix(Terms, data)
        check.nrows(expected.nrows, nrow(x), nrow(object$fitted.values), Callers.name)
        intercept <- match("(Intercept)", colnames(x), nomatch=0)
        if(intercept)
            x <- x[, -intercept, drop=FALSE]    # silently discard intercept
    }
    if(nrow(x) == 0)
        stop0("empty model matrix")
    # Fix: April 2010, allow earth to play nicely with fda with factors in x
    if(ncol(x) > ncol(object$dirs))                      # too many columns?
        x <- x[, colnames(x) %in% colnames(object$dirs), drop=FALSE] # select only the columns in dirs
    check.expanded.ncols(x, object)
    x
}
# Called by predict.earth and can also be called by users directly.
# Return object$bx if all x, subset, which.terms equal NULL.

model.matrix.earth <- function(     # returns bx
    object       = stop("no 'object' argument"),
    x            = NULL,            # x arg must not yet be expanded
    subset       = NULL,            # not used by the earth code
    which.terms  = NULL,            # not used by the earth code
    ...,                            # unused, for generic method comparibility
    env          = parent.frame(),
    trace        = 0,
    Callers.name = "model.matrix.earth") # caller's name for trace messages
{
    warn.if.dots(...)
    check.classname(object, substitute(object), "earth")
    trace <- as.numeric(check.numeric.scalar(trace, logical.ok=TRUE))
    if(is.null(x) && is.null(subset) && is.null(which.terms)) {
        if(trace >= 1)
            cat0(Callers.name, ": returning object$bx\n")
        return(object$bx)
    }
    x <- get.earth.x(object, data=x, env, trace, paste("get.earth.x from", Callers.name))
    if(is.null(which.terms))
        which.terms <- object$selected.terms
    if(!is.null(subset)) {
        # duplicates are allowed in subsets so user can specify a bootstrap sample
        check.index("subset", subset, x, allow.dups=TRUE, allow.zeros=TRUE)
        x <- x[subset, , drop=FALSE]
    }
    get.bx(x, which.terms, object$dirs, object$cuts)
}
# Called by update.earth and get.earth.x
#
# Which x should we use? The precedence is [1] the x parameter, if any,
# in this call to update [2] the $x in the earth object (which exists
# if keepxy=TRUE was used the original call to earth) [3] the x found
# in the original call to earth.
# Same applies for y, subset, weights, and wp.
# The "arg" argument is from the current call to update or predict

get.update.arg <- function(arg, argname, object, env,
                           trace1, Callers.name="update.earth", print.trace=TRUE,
                           reeval=TRUE) # TODO hack to re-evaluate
{
    if(!print.trace) # print.trace arg prevents recursion issues with trace
        trace1 = FALSE
    if(is.null(arg)) {
        temp <- try(eval(object[[argname, exact=TRUE]], envir=env), silent=TRUE)
        if(!is.null(temp) && !is.try.err(temp)) {
            if(reeval)
                arg <- object[[argname, exact=TRUE]]
            else
                arg <- temp
            if(trace1 >= 1)
                cat0(Callers.name, ": using ",
                     NROW(temp), " by ", NCOL(temp), " ", argname,
                     " saved by keepxy in original call to earth\n")
        } else {
            temp <- try(eval(object$call[[argname, exact=TRUE]], envir=env), silent=TRUE)
            if(!is.null(temp) && !is.try.err(temp)) {
                if(reeval)
                    arg <- object$call[[argname, exact=TRUE]]
                else
                    arg <- temp
                if(trace1 >= 1)
                    cat0(Callers.name, ": using ",
                         NROW(temp), " by ", NCOL(temp), " ", argname,
                         " argument from original call to earth\n")
             }
        }
    }
    arg
}
