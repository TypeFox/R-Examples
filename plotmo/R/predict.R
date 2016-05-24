# predict.R: plotmo wrapper functions for predict()

# Returns an n x 1 matrix (unless nresponse=NULL then returns an n x q dataframe)
#
# The newdata argument can be a positive integer n, which is the same as
# newdata=NULL but may return only n rows if that is more efficient.
# This is for efficiency in plotmo_meta.

plotmo_predict <- function(object, newdata, nresponse,
                           type, expected.levs, trace, inverse.func=NULL, ...)
{
    # handle special case where newdata specifies the number of rows
    nrows <- 0
    if(is.numeric(newdata) && length(newdata) == 1 && newdata > 0) {
        nrows <- newdata
        newdata <- NULL
    }
    if(is.null(newdata)) {
        # It generally faster to use newdata=NULL. But not all models process
        # the type argument with null newdata. So here we check for some models
        # that are known good that way.  The inherits function is not used here
        # because for example a glm model inherits("lm") but with NULL newdata
        # doesn't process type as we might hope.
        if(class(object)[1] %in% c("lm", "earth"))
            trace2(trace, "calling predict.%s with NULL newdata\n", class(object)[1])
        else { # assume object cannot handle newdata=NULL
            trace2(trace, "plotmo_predict with NULL newdata%s, %s",
                   if(nrows) sprintf(" (nrows=%d)", nrows) else "",
                   "using plotmo_x to get the data\n")
            newdata <- plotmo_x(object, trace)
            if(nrows)
                newdata <- newdata[seq_len(nrows),,drop=FALSE]
            trace2(trace,
                   "will use the above data instead of newdata=NULL for predict.%s\n",
                   class(object)[1])
        }
    } else
        print_summary(newdata, "newdata", trace)

    yhat <- plotmo.predict(object=object, newdata=newdata, type=type, ...,
                           TRACE=if(trace >= 2) trace else trace.call.global)

    temp <- process.y(yhat, object, type, nresponse,
                      expected.len=nrow(newdata), expected.levs, trace,
                      fname="predict")

    yhat <- apply.inverse.func(inverse.func, temp$y, object, trace)

    list(yhat       = yhat, # n x 1 matrix (unless nresponse=NULL then an n x q dataframe)
         resp.levs  = temp$resp.levs,
         resp.class = temp$resp.class)
}
# TRACE is passed to do.call.trace (if TRACE>0 print the call to predict)
plotmo.predict <- function(object, newdata, type, ..., TRACE)
{
    stopifnot.string(type)
    UseMethod("plotmo.predict")
}
# this handles a common mistake
plotmo.predict.list <- function(object, newdata, type, ..., TRACE)
{
    stop0("object does not have a predict method")
}
# plotmo.predict.default calls predict for the given object,
# and does tracing and error handling.
#
# It also allows use to pre-program default args for predict,
# which can be overruled or augmented by args passed in dots.
# These defaults args must be specified in the calling function.  For example
#   plotmo.predict.default(object, newdata, type=type, def.foo=3, ...)
# will pass foo=3 to predict --- unless the caller of plotmo passes
# predict.foo=0 to plotmo, which will override the default and pass foo=0
# to predict.
# When specifying defaults, use the full arg name (no abbreviations)
# prefixed by "def.".

plotmo.predict.default <- function(
    object,
    newdata,
    ...,       # extra args to predict, first typically is type="xxx"
    TRACE,     # passed to do.call.trace (if TRACE>0 print the call to predict)
    FUNC=NULL) # predict function, NULL means use stats::predict
{
    fname <- "PREDICTFUNC"
    if(is.null(FUNC)) {
        FUNC <- stats::predict
        fname <- "stats::predict"
    }
    # Create arg list for predict.
    # We invoke deprefix directly (and not call.dots) because we have to
    # specify a DROP argument and also do a bit of other processing.
    # OBJECT and NEWDATA must be passed as unnamed arguments to predict,
    # because different predict methods use different arg names for these.
    # We want to allow the user to pass normal (unprefixed) dots argument to
    # predict.  So here we use KEEP=NULL but drop any plot arguments, and
    # any prefixed dot arguments that are necessary elsewhere in plotmo.
    # We can't specify a FUNC argument to deprefix because we don't
    # know which specific predict.method will be called (a few lines down).

    args <- deprefix(FUNC=NULL,
                DROP=paste0("w1. SHOWCALL FORCEPREDICT PLOT.ARGS PAR.ARGS PLOTMO.ARGS"),
                PREFIX="predict.", FNAME=fname,
                force.anon1=object, force.anon2=newdata, ...)

    yhat <- do.call.trace(func=FUNC, args=args, fname=fname, trace=TRACE)

    if(is.null(yhat) || length(yhat) == 0)
        stopf("failed call to predict(%s)", list.as.char(args))

    yhat  # plausibility of yhat will be checked shortly in plotmo_predict
}
# Like plotmo.predict.default but first convert newdata to a matrix.
# Needed because some predict methods require a matrix, not a data.frame.

plotmo.predict.defaultm <- function(object, newdata, type, ..., TRACE, FUNC=NULL)
{
    stopifnot(is.data.frame(newdata))
    check.df.numeric.or.logical(newdata)
    # following calls predict.xxx where xxx is the class of object
    plotmo.predict.default(object, data.matrix(newdata), type=type, ..., TRACE=TRACE)
}
