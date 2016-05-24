# lib.R: miscellaneous functions for both earth and plotmo
#        functions in this file are in alphabetical order

any1 <- function(x)
{
    any(x != 0) # like any but no warning if x not logical
}
cat0 <- function(...) # cat with no added spaces
{
    cat(..., sep="")
}
check <- function(object, object.name, check.name, check.func, na.ok=FALSE)
{
    any <- check.func(object)
    if(na.ok)
        any <- any[!is.na(any)]
    else {
        which.na <- which(is.na(any))
        if(length(which.na)) {
            stopf("NA in %s\n       %s[%d] is %g",
                  object.name, unquote(object.name),
                  which.na[1], object[which.na[1]])
        }
    }
    if(any(any)) {
        which <- which(check.func(object))
        stopifnot(length(which) > 0)
        stopf("%s in %s\n       %s[%d] is %g",
              check.name, object.name, unquote(object.name),
              which[1], object[which[1]])
    }
}
# TODO commented out the following because it is too slow for big data
#      (the as.character is very slow)
#
# # The args argument is assumed to be a list of arguments for do.call.
# # An argument in args will be an unforced promise if it couldn't be
# # evaluated earlier e.g. if call.plot was invoked with arg=nonesuch.
# # If an argument is such an unforced promise, issue an error message now
# # to prevent very confusing error messages later.  To do this, we have to
# # determine if the arg is a promise, which we do with the if statement
# # below.
# # This makes me nervous, because the R language manual says "There is
# # generally no way in R code to check whether an object is a promise or not".
#
# check.do.call.args <- function(func, args, fname)
# {
#     stopifnot(is.list(args))
#     for(i in seq_along(args)) {
#         if(length(args[i]) == 1 && !is.na(args[i]) &&
#            substr(as.character(args[i]), 1, 2) == "..") {
#             printf("\n")
#             s <- paste0(strwrap(list.as.char(args),
#                         width=getOption("width"), exdent=7), collapse="\n")
#             stop0("cannot evaluate '", names(args)[i],
#                   "' in\n       ", fname, "(", s, ")")
#         }
#     }
# }

# mostly for checking user arguments (so error wording is for that)
# but also occasionally used for other sanity checking

check.boolean <- function(b) # b==0 or b==1 is also ok
{
    if(length(b) != 1)
        stop0("the ", short.deparse(substitute(b), "given"),
              " argument is not FALSE, TRUE, 0, or 1")
    if(!(is.logical(b) || is.numeric(b)) || is.na(b) || !(b == 0 || b == 1))
        stop0(short.deparse(substitute(b), "the argument"), "=", b,
            " but it should be FALSE, TRUE, 0, or 1")
    b != 0 # convert to logical
}
check.classname <- function(object, substituted.object, expected.classname)
{
    if(is.null(object))
        stopf("object is NULL (expected an object of class \"%s\")",
              expected.classname)
    if(!inherits(object, expected.classname)) {
        stopf("the class of '%s' is \"%s\" but expected the class to be \"%s\"",
              paste.trunc(substituted.object, maxlen=30),
              class(object)[1], expected.classname)
    }
}
check.integer.scalar <- function(object, min=NULL, max=NULL, null.ok=FALSE,
                                 na.ok=FALSE, logical.ok=TRUE,
                                 char.ok=FALSE,
                                 object.name=short.deparse(substitute(object)))
{
    stop.msg <- function()
    {
        s.null    <- if(null.ok)    ", or NULL"          else ""
        s.na      <- if(na.ok)      ", or NA"            else ""
        s.logical <- if(logical.ok) ", or TRUE or FALSE" else ""
        s.char    <- if(char.ok)    ", or a string"      else ""
        stop0(object.name, "=", object[1], " but it should be an an integer",
              s.null, s.na, s.logical, s.char)
    }
    if(is.character(object)) {
        if(!char.ok || length(object) != 1)
            stop.msg()
    } else {
        check.numeric.scalar(object, null.ok, na.ok, logical.ok,
                             char.ok.msg=char.ok, object.name=object.name)
        if(!is.null(object) && !is.na(object)) {
            if(object != floor(object))
                stop.msg()
            if(!is.null(min) && !is.null(max) && (object < min || object > max)) {
                stop0(object.name, "=", object,
                      " but it should be between ", min, " and ", max)
            }
            if(!is.null(min) && object < min) {
                stop0(object.name, "=", object,
                      " but it should be at least ", min)
            }
            if(!is.null(max) && object > max) {
                stop0(object.name, "=", object,
                      " but it should not be greater than ", max)
            }
        }
    }
    object
}
check.level.arg <- function(level, zero.ok)
{
    if(anyNA(level) || is.null(level)) # treat NA and NULL as 0
        level <- 0
    check.numeric.scalar(level)
    if(!((zero.ok && level == 0) || level >= .5 || level < 1)) {
        stop0("level=", level, " but it should be ",
              if(zero.ok) "zero or " else "", "between 0.5 and 1")
    }
    level
}
check.no.na.in.mat <- function(object)
{
    if(anyNA(object)) { # quick initial check
        # detailed check for detailed error message
        for(icol in seq_along(ncol(object))) {
            check.name <-
                if(!is.null(colnames(object)))
                    colnames(object)[icol]
                else
                    sprintf("%s[,%d]",
                        short.deparse(substitute(object), "matrix"), icol)

            check(object, check.name, "NA", is.na, na.ok=FALSE)
        }
    }
}
# x can be a data.frame or matrix
check.df.numeric.or.logical <- function(x, xname=trunc.deparse(substitute(x)))
{
    stopifnot(!is.null(x), length(dim(x)) == 2)
    for(icol in seq_len(NCOL(x))) {
        if(!is.numeric(x[,icol]) && !is.logical(x[,icol]))
            stopf("the class of %s is \"%s\" (expected numeric or logical)",
                  colname(x, icol, xname), class(x[,icol]))
        is.na <- is.na(x[,icol])
        if(any(is.na))
            stopf("%s[%g] is NA", colname(x, icol, xname), which(is.na)[1])
        is.infinite <- !is.finite(x[,icol])
        if(any(is.infinite))
            stopf("%s[%g] is Inf", colname(x, icol, xname), which(is.infinite)[1])
    }
}
check.numeric.scalar <- function(object, null.ok=FALSE,
                                 na.ok=FALSE, logical.ok=FALSE,
                                 char.ok.msg=FALSE, # only affects error msg
                                 object.name=short.deparse(substitute(object)))
{
    s.logical <- if(logical.ok) ", or TRUE or FALSE" else ""
    if(na.ok)
        logical.ok <- TRUE # needed because NA is a logical
    any.na <- !is.null(object) && anyNA(object)
    if(is.null(object)) {
        if(!null.ok)
            stop0(object.name, "=NULL is not allowed")
    } else if(any.na && !na.ok)
        stop0(object.name, "=NA is not allowed")
    else if(!is.numeric(object) && !(is.logical(object) && logical.ok)) {
        s.na   <- if(na.ok)       ", or NA"       else ""
        s.null <- if(null.ok)     ", or NULL"     else ""
        s.char <- if(char.ok.msg) ", or a string" else ""
        stopf("'%s' must be numeric%s%s%s%s (whereas its current class is \"%s\")",
              object.name, s.null, s.na, s.char, s.logical, class(object)[1])
    }
    else if(length(object) != 1) {
        stopf("the length of '%s' must be 1 (whereas its current length is %d)",
              object.name, length(object))
    }
    object
}
# We allow 20% of x to be nonpositive, useful if the response is essentially
# positive, but the predicted response has a few nonpositive values at the extremes.
# Needed for example if we will later take log(x) or sqrt(x).

check.that.most.are.positive <- function(x, xname, user.arg, non.positive.msg, frac.allowed=.2)
{
    check.numeric.scalar(frac.allowed)
    stopifnot(frac.allowed >= 0, frac.allowed <= 1)
    nonpos <- x <= 0
    if(sum(nonpos, na.rm=TRUE) > frac.allowed * length(x)) { # more than frac.allowed nonpos?
        ifirst <- which(nonpos)[1]
        stop0(sprintf(
                "%s is not allowed because too many %ss are %s\n",
                user.arg, unquote(xname), non.positive.msg),
              sprintf(
                "       %.2g%% are %s (%g%% is allowed)\n",
                100 * sum(nonpos) / length(x), non.positive.msg, 100 * frac.allowed),
               sprintf("       e.g. %s[%d] is %g", unquote(xname), ifirst, x[ifirst]))
    }
}
check.vec <- function(object, object.name, expected.len=NA, logical.ok=TRUE, na.ok=FALSE)
{
    if(!(NROW(object) == 1 || NCOL(object) == 1))
        stop0(object.name, " is not a vector\n       ",
              object.name, " has dimensions ", NROW(object), " by ", NCOL(object))
    if(!((logical.ok && is.logical(object)) || is.numeric(object)))
        stop0(object.name, " is not numeric")
    if(!is.na(expected.len) && length(object) != expected.len)
        stop0(object.name, " has the wrong length ",
              length(object), ", expected ", expected.len)
    if(na.ok)
        object[is.na(object)] <- 1 # prevent check is.finite from complaining
    else
        check(object, object.name, "NA", is.na)
    check(object, object.name, "non-finite value", function(object) {!is.finite(object)})
}
# returns the column name, if that is not possible then something like x[,1]
colname <- function(object, i, object.name=trunc.deparse(substitute(object)))
{
    check.numeric.scalar(i)
    check.index(i, object.name, object, is.col.index=TRUE, allow.negatives=FALSE)
    colnames <- safe.colnames(object)
    if(!is.null(colnames))
        colnames[i]
    else if(NCOL(object) > 1)
        sprintf("%s[,%g]", object.name, i)
    else
        sprintf(object.name)
}
# if trace>0 or the func fails, then print the call to func

do.call.trace <- function(func, args,
                          fname=short.deparse(deparse(func), "FUNC"), trace=0)
{
    stopifnot(is.logical(trace) || is.numeric(trace), length(trace) == 1)
    # TODO commented out the following because it is too slow for big data
    # check.do.call.args(func, args, fname)
    trace <- as.numeric(trace)
    if(trace > 0)
        printf.wrap("%s(%s)\n", fname, list.as.char(args))
    try <- try(do.call(what=func, args=args), silent=TRUE)
    if(is.try.err(try)) {
        if(trace == 0) # didn't print call above? then print it now
            printf.wrap("\n%s(%s)\n\n", fname, list.as.char(args))
        else if(trace >= 2) # TODO is this best?
            printf("\n")
        # Re-call func so user can do a traceback within the function.  Note that
        # if do.call.trace was called with try, this will be caught by that try.
        # TODO is there a better way to achieve this, perhaps using tryCatch
        #      this could be confusing if func has side effects (unlikely)
        do.call(what=func, args=args)
        # should never get here
        stop0("second do.call(", fname,
              ", ...) did not give the expected error: ", try[1])
    }
    invisible(try) # TODO is invisible necessary?
}
# identical to base::eval() but has trace and expr.name arguments
eval.trace <- function(
    expr,
    envir     = parent.frame(),
    enclos    = if(is.list(envir) || is.pairlist(envir))
                   parent.frame()
                else
                   baseenv(),
    trace     = 0,
    expr.name = NULL)
{
    stopifnot(is.environment(envir))
    stopifnot(is.environment(enclos))
    if(trace >= 2)
        printf("eval(%s, %s)\n",
            if(is.null(expr.name)) trunc.deparse(substitute(expr))
            else                   expr.name,
            environment.as.char(envir))

    eval(expr, envir, enclos)
}
# This function is used for xlim and ylim.
# If xlim[1] == xlim[2], then plot() issues a message.  We don't want that,
# so use this function to make sure xlim[2] is different to xlim[1].

fix.lim <- function(lim)
{
    if(!is.null(lim)) {
        stopifnot(is.numeric(lim), length(lim) == 2)
        # constants below are arbitrary
        small <- max(1e-6, .001 * abs(lim[1]), .001 * abs(lim[2]))
        if(abs(lim[2] - lim[1]) < small) # illegal lim?
            lim <- c(lim[1] - small, lim[2] + small)
    }
    lim
}
# Ensure all columns of x have column names.  Won't overwrite existing
# column names (TODO except possibly if make.unique kicks in).

gen.colnames <- function(x, prefix="x", alt.prefix=prefix, trace=0, xname=NULL)
{
    if(NCOL(x) == 0)
        return(NULL)
    # If prefix is long and has characters like ( or [ then use the
    # alternate prefix.  This is sometimes necessary when prefix is
    # generated using deparse and the arg is something like
    # "cbind(trees$Volume,trees$Volume+100)"
    if(any(nchar(prefix) > 30) && grepany("[([,]", prefix)) {
        trace2(trace, "using alt.prefix \"%s\" instead of prefix \"%s\"\n",
               alt.prefix, prefix)
        prefix <- alt.prefix
    }
    stopifnot(length(prefix) <= NCOL(x))
    prefix <- substr(prefix, 1, 60)
    new.colnames <-
        if(NCOL(x) == length(prefix))
            prefix
        else if(grepany("\\[", prefix))
            new.colnames <- paste0(prefix, "[", seq_len(NCOL(x)), "]")
        else
            new.colnames <- paste0(prefix, seq_len(NCOL(x)))
    colnames <- org.colnames <- colnames(x)
    if(is.null(colnames))
        colnames <- new.colnames
    else {
        missing <- !nzchar(colnames)
        if(any(missing))
            colnames[missing] <- new.colnames[missing]
    }
    colnames <- make.unique(strip.space(colnames))
    if(trace >= 2 && !identical(org.colnames, colnames))
        trace2(trace, "%s colname%s %s now %s\n",
            if(is.null(xname)) trunc.deparse(substitute(x)) else xname,
            if(length(colnames) > 1) "s were" else " was",
            if(is.null(org.colnames)) "NULL"
            else paste.trunc(quotify(org.colnames)),
            paste.trunc(quotify(colnames)))
    colnames
}
get.mean.rsq <- function(rss, tss, wp)
{
    if(is.null(wp))
        wp <- repl(1, length(rss))
    stopifnot(length(rss) == length(tss), length(wp) == length(tss))
    total.rsq <- 0
    for(iresp in seq_along(rss))
        total.rsq <- total.rsq + wp[iresp] * get.rsq(rss[iresp], tss[iresp])
    sum(total.rsq) / sum(wp)
}
# Get the environment for evaluating the model data:
#   1. Return the environment in which the model function
#      was originally called.
#   2. Else if the model already has an attribute .Environment, use that.
#   3. Else return the environment in which the caller of this function
#      was called (e.g. return the environment of plotmo's caller).

get.model.env <- function(object, object.name="object", trace=0)
{
    # check args, because this func is called very early in plotmo (and friends)
    stopifnot.string(object.name)
    check.numeric.scalar(trace, logical.ok=TRUE)
    if(!is.list(object))
        stopf("%s is not an S3 model", object.name)
    if(class(object)[1] == "list")
        stopf("%s is a plain list, not an S3 model", object.name)

    # following will fail for non-formula models because they have no terms field
    terms <- try(terms(object), silent=trace < 3)
    if(!is.try.err(terms) && !is.null(terms)) {
        model.env <- attr(terms, ".Environment")
        if(is.null(model.env))
            stop0("attr(terms, \".Environment\") is NULL")
        if(!is.environment(model.env))
            stop0("attr(terms, \".Environment\") is not an environment")
        else {
            trace2(trace, "using the environment saved with the %s model: %s\n",
                   class(object)[1], environment.as.char(model.env))
            return(model.env)
        }
    }
    model.env <- attr(object, ".Environment")
    if(is.environment(model.env)) {
        trace2(trace, "using attr(object,\".Environment\") saved with %s model: %s\n",
               class(object)[1], environment.as.char(model.env))
        return(model.env)
    }
    if(!is.null(model.env))
        stop0("attr(object, \".Environment\") is not an environment")

    model.env <- parent.frame(n=2) # caller of the function that called model.env
    if(trace >= 2) {
        callers.name <- callers.name()
        printf("assuming the environment of the %s model is that of %s's caller: %s\n",
               class(object)[1], callers.name, environment.as.char(model.env))
    }
    return(model.env)
}
get.rsq <- function(rss, tss)
{
    rsq <- 1 - rss / tss
    # following makes testing easier across machines in presence of numerical error
    rsq[rsq > -1e-5 & rsq < 1e-5] <- 0
    rsq
}
get.weighted.rsq <- function(y, yhat, w=NULL) # NAs will be dropped before calc
{
    stopifnot(length(y) > 0, length(y) == length(yhat))
    if(is.null(w)) {
        is.na <- is.na(y) | is.na(yhat)
        y    <- y[!is.na]
        yhat <- yhat[!is.na]
        if(length(y) == 0)
            stop0("length(y) == 0 after deleting NAs in y or yhat")
        rss <- sum((y - yhat)^2)
        tss <- sum((y - mean(y))^2)
    } else {
        stopifnot(length(w) == length(yhat))
        is.na <- is.na(y) | is.na(yhat) | is.na(w)
        y    <- y[!is.na]
        yhat <- yhat[!is.na]
        w    <- w[!is.na]
        if(length(y) == 0)
            stop0("length(y) == 0 after deleting NAs in y or yhat or w")
        rss <- sum(w * (y - yhat)^2)
        tss <- sum(w * (y - weighted.mean(y, w))^2)
    }
    get.rsq(rss, tss)
}
# TRUE if pattern is in any of the strings in x
grepany <- function(pattern, x, ignore.case=FALSE, ...)
{
    any(grepl(pattern, x, ignore.case=ignore.case, ...))
}
# returns an index, choices is a vector of strings
imatch.choices <- function(arg, choices,
            argname=short.deparse(substitute(arg), "function"),
            err.msg.has.index=FALSE) # TRUE if integer "arg" is legal elsewhere
{
    if(!is.character(arg) || length(arg) != 1 || !nzchar(arg))
         stopf("illegal '%s' argument\nChoose%s one of: %s",
               argname,
               if(err.msg.has.index) " an integer index or" else "",
               quotify(choices))
    imatch <- pmatch(arg, choices)
    if(anyNA(imatch)) {
        imatch <- NULL
        for(i in seq_along(choices))
            if(pmatch(arg, choices[i], nomatch=0))
                imatch <- c(i, imatch)
        if(length(imatch) == 0) {
            if(length(choices) == 1)
                stopf("%s=\"%s\" is not allowed\n       Only%s %s is allowed",
                      argname, paste(arg),
                      if(err.msg.has.index) " 1 or" else "",
                      quotify(choices))
            else
                stopf("%s=\"%s\" is not allowed\nChoose%s one of %s",
                      argname, paste(arg),
                      if(err.msg.has.index) " an integer index or" else "",
                      quotify(choices))
        }
        if(length(imatch) > 1)
            stopf("%s=\"%s\" is ambiguous\nChoose%s one of %s",
                  argname, paste(arg),
                  if(err.msg.has.index) " an integer index or" else "",
                  quotify(choices))
    }
    imatch
}
is.integral <- function(object)
{
    all(floor(object) == object)
}
# is.specified's main purpose is to see if a plot component should be
# drawn, i.e., to see if the component "has a color"

is.specified <- function(object)
{
    !is.null(object) && !anyNA(object) && !identical(object, 0)
}
is.try.err <- function(object)
{
    class(object)[1] == "try-error"
}
# returns the expanded arg
match.arg1 <- function(arg, argname=deparse(substitute(arg)))
{
    formal.args <- formals(sys.function(sys.parent()))
    formal.argnames <- eval(formal.args[[argname]])
    formal.argnames[imatch.choices(arg[1], formal.argnames, argname)]
}
# returns a string, choices is a vector of strings
match.choices <- function(arg, choices, argname=deparse(substitute(arg)))
{
    choices[imatch.choices(arg, choices, argname)]
}
# This uses the object's .Environment attribute, which was
# pre-assigned to the object via get.model.env
# If this gives an error saying that class(model.env) is "NULL"
# then that pre-assignment wasn't done.

model.env <- function(object)
{
    model.env <- attr(object, ".Environment")
    if(!is.environment(model.env))
        stopf("class(model.env) is \"%s\"", class(model.env)[1])
    model.env
}
# Like as.data.frame() but retains the original colnames, if any, and can
# handle matrices from the Matrix etc. packages, if as.matrix() works for
# them.  Also it has a stringsAsFactors argument which works even if x is
# already a data.frame.

my.data.frame <- function(x, trace, stringsAsFactors=TRUE)
{
    if(is.data.frame(x)) {
        if(stringsAsFactors) {
            # Convert any character columns to factors.  Note as.data.frame
            # won't do this for us when x is already a data.frame.
            # We don't have a levels argument to pass to factor()
            # but I believe that this will not be a problem in the
            # context in which we use my.data.frame (plotmo_x).
            for(i in seq_len(length(x)))
                if(is.character(x[[i]]))
                    x[[i]] <- factor(x[[i]])
            }
       return(x)
    }
    df <- try(as.data.frame(x, stringsAsFactors=stringsAsFactors), silent=TRUE)
    if(is.try.err(df)) {
        # come here for sparse matrices from the Matrix package
        df <- try(as.matrix(x))
        if(is.try.err(df))
            stopf("Could not convert '%s' object to a data.frame or matrix",
                  class(x)[1])
        df <- as.data.frame(df, stringsAsFactors=stringsAsFactors)
        trace2(trace, "converted %s object to data.frame\n", class(x)[1])
    }
    colnames(df) <- safe.colnames(x) # restore original column names
    df
}
my.fixed.point <- function(x, digits)
{
    if(NROW(x) > 2) # only use fixed point if more than intercept and one other term
        x <- apply(x, 2, zapsmall, digits+1)
    x
}
# If s is a string vector s, return the number of lines in
# the element that has the most lines
# Examples: nlines(c(" ", " \n ") is 2
#           nlines(c(" ", " \n")  is 2
#           nlines(" ")           is 1
#           nlines("")            is 0  (special case)

nlines <- function(s)
{
    if(!nzchar(s[1])) # special case, caption="" is not printed
        0
    else
        length(strsplit(s, "\n")[[1]])
}
paste.c <- function(object, maxlen=16) # return "x1" or "c(x1, x2)"
{
    if(length(object) == 1)
        paste.trunc(object)
    else
        paste0("c(", paste.trunc(object, collapse=",", maxlen=maxlen), ")")
}
paste.collapse <- function(...)
{
    paste(..., collapse=" ")
}
# collapse, and truncate if strings in ... are too long
paste.trunc <- function(..., sep=" ", collapse=" ", maxlen=60)
{
    s <- paste(..., sep=sep, collapse=collapse)
    if(nchar(s) > maxlen) {
        stopifnot(maxlen > 3)
        s <- paste0(substr(s, 1, maxlen-3),
                    if(substr(s, maxlen-3, maxlen-3) == ".") ".." # avoid 4 dots
                    else                                     "...")
    }
    s
}
pastef <- function(s, fmt, ...) # paste the printf style args to s
{
    paste0(s, sprintf(fmt, ...))
}
pt.cex <- function(ncases, npoints=ncases)
{
    n <- if(npoints > 0) min(npoints, ncases) else ncases

    if     (n >= 20000) .2
    else if(n >=  5000) .3
    else if(n >=  3000) .4
    else if(n >=  1000) .6
    else if(n >=  300)  .8
    else if(n >=  30)   1
    else                1.2
}
print.first.few.elements.of.vector <- function(x, trace, name=NULL)
{
    try(cat(" min", min(x), "max", max(x)), silent=TRUE)
    spaces <- "               "
    if(!is.null(name))
        spaces <- sprintf("%*s", nchar(name), " ")  # nchar spaces
    cat0("\n", spaces, " value")
    len <- if(trace >= 4) length(x)
           else           min(if(is.logical(x)) 20 else 10, length(x))
    if(is.logical(x))
        for(i in 1:len)
            cat0(if(x[i]) " T" else " F")
    else
        for(i in 1:len)
            cat0(" ", x[i])
    if(length(x) > len)
        cat(" ...")
   cat("\n")
    if(trace >= 4) {
        cat("\n")
        print(summary(x))
    }
}
printf <- function(fmt, ...) # like c printf
{
    cat(sprintf(fmt, ...), sep="")
}
# like printf but wrap at terminal width
# exdent=NULL for automatic determination of xdent (line up to func opening paren)
# TODO maxlen seems to be ignored, strwrap truncates before that?
printf.wrap <- function(fmt, ..., exdent=NULL, maxlen=2000)
{
    s <- paste.trunc(paste.collapse(sprintf(fmt, ...)), maxlen=maxlen)
    if(is.null(exdent)) {
        # align to opening paren of func call e.g. "graphics::par(xxx)" or "foo$method("
        # TODO this doesn't account for leading newlines if any
        exdent <- 4
        igrep <- gregexpr("[ ._$:[:alnum:]]+\\(", s)[[1]]
        if(igrep[1] == 1) {
            len <- attr(igrep, "match.length")[1]
            exdent <- min(25, len)
        }
    }
    # strwrap doesn't preserve newlines in the input string, so do it manually :(
    for(i in seq_len(nchar(s))) # print leading newlines
        if(substr(s, i, i) == "\n") cat0("\n") else break

    cat(paste0(strwrap(s, width=getOption("width"), exdent=exdent),
               collapse="\n"))

    if(nchar(s) > i) for(j in nchar(s):i) # print trailing newlines
        if(substr(s, j, j) == "\n") cat0("\n") else break
}
# like short.deparse but quotify the deparsed obj (unless the alternative is used)
quote.deparse <- function(object, alternative="object")
{
    s <- strip.deparse(object)
    if(nchar(s) > 60)
        alternative
    else
        quotify(s, quote="'")
}
quote.with.c <- function(names) # return "x" or c("x1", "x2")
{
    if(length(names) == 1)
        sprintf("\"%s\"", names)
    else
        sprintf("c(%s)", paste0("\"", paste(names, collapse="\", \""), "\""))
}
quotify <- function(s, quote="\"") # add quotes and collapse to a single string
{                                  # called quotify because quote is taken
    if(is.null(s))
        "NULL"
    else if(length(s) == 0)
        paste0(quote, quote)       # not sure what is best here
    else if(substr(s[1], 1, 1) == quote) # already has quotes?
        paste.collapse(s)
    else
        paste0(quote, paste(s, collapse=paste0(quote, " ", quote)), quote)
}
# like quotify, but use the alternative name if s is too long
quotify.short <- function(s, alternative="object", quote="\"")
{
    stopifnot(is.character(s))
    s <- paste0(s, collapse="")
    if(nchar(s) > 60) # 60 is arb but seems ok for plot titles etc
        alternative
    else
        quotify(s, quote)
}
quotify.trunc <- function(s, quote="\"", maxlen=60)
{
    stopifnot(is.character(s))
    s <- quotify(s, quote)
    if(nchar(s) > maxlen) {
        stopifnot(maxlen > 3)
        paste0(substr(s, 1, maxlen-3), "...")
    } else
        s
}
range1 <- function(object, ...)
{
    stopifnot(length(dim(object)) <= 2)
    if(!is.null(dim(object)))
        object <- object[,1]
    if(is.factor(object))
        c(1, nlevels(object))
    else
        range(object, finite=TRUE, ...)
}
repl <- function(object, length.out)
{
    check.numeric.scalar(length.out)
    stopifnot(floor(length.out) == length.out)
    stopifnot(length.out > 0)
    rep(object, length.out=length.out)
}
# the standard colnames() can crash for certain objects
# TODO figure out when and why

safe.colnames <- function(object)
{
    colnames <- try(colnames(object), silent=TRUE)
    if(is.try.err(colnames))
        NULL
    else
        colnames
}
# if deparse(object) is too long, return the alternative
short.deparse <- function(object, alternative="object")
{
    s <- strip.deparse(object)
    if(nchar(s) > 60)
        alternative
    else
        s
}
stop0 <- function(...)
{
    stop(..., call.=FALSE)
}
stopf <- function(fmt, ...) # args like printf
{
    stop(sprintf(fmt, ...), call.=FALSE)
}
# stop if s is not a one element character vector
stopifnot.string <- function(s, name=short.deparse(substitute(s)),
                             null.ok=FALSE, allow.empty=FALSE)
{
    if(is.null(s)) {
        if(null.ok)
            return()
        else
            stop0(name, " is NULL (it should be a string)")
    }
    if(!is.character(s))
        stop0(name, " is not a character variable (class(",
              name, ") is \"", class(s), "\")")
    if(length(s) != 1)
        stop0(name, " has more than one element\n       ",
              name, " = c(", paste.trunc("\"", s, "\"", sep=""), ")")
    if(!allow.empty && !nzchar(s))
        stop0(name, " is an empty string")
}
strip.deparse <- function(object) # deparse, collapse, remove most white space
{
    s <- strip.space(paste0(deparse(object), collapse=""))
    gsub(",", ", ", s) # put back space after commas
}
strip.space <- function(s)
{
    gsub("[ \t\n]", "", s)
}
# like text, but with a white background
# TODO sign of adj is backwards?

text.on.white <- function(x, y, label,
                          cex=1, adj=.5, font=1, xmar=.3, srt=0, white="white", ...)
{
    stopifnot(length(label) == 1)
    if(length(adj) == 1)
        adj <- c(adj, .5)
    width       <- strwidth(label,  cex=cex, font=font)
    char.width  <- strwidth("X",    cex=cex, font=font)
    height      <- strheight(label, cex=cex, font=font)
    char.height <- strheight("X",   cex=cex, font=font)
    if(srt == 0) {
        rect(x - adj[1]     * width  - xmar * char.width,
             y - adj[2]     * height - .3   * char.height, # .3 for extra space at bottom
             x + (1-adj[1]) * width  + xmar * char.width,
             y + (1-adj[2]) * height + .1   * char.height,
             col=white, border=NA)
        text(x=x, y=y, labels=label, cex=cex, adj=adj, font=font, ...)
    } else if(srt == 90 || srt == -90) {
        # width and height are in usr coords, adjust these for flip of coords
        usr <- par("usr") # xmin, xmax, ymin, ymax
        xrange <- abs(usr[2] - usr[1])
        yrange <- abs(usr[4] - usr[3])
        height <- xrange / yrange * height
        width  <- yrange / xrange * width
        char.height <- xrange / yrange * char.height
        char.width  <- yrange / xrange * char.width
        rect(x + (1-adj[1]) * height,                    # left
             y + (1-adj[2]) * width + xmar * char.width, # bottom
             x - adj[1]     * height,                    # right
             y - adj[2]     * width - xmar * char.width, # top
             col=white, border=NA)
        text(x=x, y=y, labels=label, cex=cex, adj=adj, font=font, srt=srt, ...)
    }
    else
        stop0("srt=", srt, " is not allowed (only 0, 90, and -90 are supported)")
}
to.logical <- function(object, len)
{
    xlogical <- repl(FALSE, len)
    xlogical[object] <- TRUE
    xlogical
}
trace1 <- function(trace, fmt, ...)
{
    stopifnot(is.numeric(trace))
    if(trace >= 1)
        cat(sprintf(fmt, ...), sep="")
}
trace2 <- function(trace, fmt, ...)
{
    stopifnot(is.numeric(trace))
    if(trace >= 2)
        cat(sprintf(fmt, ...), sep="")
}
# Truncate deparse(object) if it is too long.
# Necessary because deparse(substitute(x)) might return something very
# long, like c(1000, 1001, 1002, 1003, 1004, 1005, 1006, 1008, 1009, etc.)
# Return a one element character vector.

trunc.deparse <- function(object, maxlen=60)
{
    s <- strip.deparse(object)
    if(nchar(s) > maxlen) {
        stopifnot(maxlen > 3)
        paste0(substr(s, 1, maxlen-3), "...")
    } else
        s
}
# Return the number of lines in s (where lines are separated by \n).
try.eval <- function(
    expr,
    envir     = parent.frame(),
    trace     = 0,
    expr.name = NULL,
    silent    = trace < 2)
{
    if(trace && is.null(expr.name))
        expr.name <- trunc.deparse(substitute(expr))
     try(eval.trace(expr, envir, trace=trace, expr.name=expr.name), silent=silent)
}
unquote <- function(s) # remove leading and trailing quotes, if any
{
    if(is.character(s)) {
        s <- gsub("^\"|^'", "", s)  # leading quotes
        s <- gsub("\"$|'$", "", s)  # trailing quotes
    }
    s
}
# warn.if.not.all.finite helps preempt confusing message from code later.
# Return TRUE if warning issued.

warn.if.not.all.finite <- function(object, text="unknown")
{
    is.factors <- sapply(object, is.factor)
    if(any(is.factors)) {
        if(NCOL(object) == 1 || all(is.factors)) # TODO suspect
            return(FALSE)
        object <- object[, !is.factors]  # remove factor columns before is.finite check
    }
    if(any(sapply(object, is.na))) {
        warning0("NA in ", text)
        return(TRUE)
    }
    if(!all(sapply(object, is.finite))) {
        warning0("non finite value in ", text)
        return(TRUE)
    }
    FALSE
}
warnf <- function(fmt, ...) # args like printf
{
    warning(sprintf(fmt, ...), call.=FALSE)
}
warning0 <- function(...)
{
    warning(..., call.=FALSE)
}
