# meta.R: plotmo function to get the "metadata" from the model

plotmo_type <- function(object, trace, fname="plotmo", type, ...)
{
    if(is.null(type)) # get default type for this object class?
        type <- plotmo.type(object, ...)
    else {
        stopifnot.string(type)
        if(pmatch(type, "terms", nomatch=0))
            stop0("type=\"terms\" is not allowed by ", fname)
    }
    type
}
plotmo_residtype <- function(object, trace, fname="plotmo", type, ...)
{
    if(is.null(type)) # get default type for this object class?
        type <- plotmo.residtype(object, ...)
    else
        stopifnot.string(type)
    type
}
# In plotmo and plotres there is some general data we need about the
# model.  For example, the response name. This routine provides that
# data, which we call "metadata".
#
# Also, plotmo and plotres should work automatically, as much as possible,
# without requiring the user to specify arguments.  This routine
# facilitates that.
#
# For example, it converts the default nresponse=NA to a sensible column
# number in the response. It will issue an error message if it can't do
# that.
#
# It also converts the default type=NULL into an appropriate
# model-specific type for predict().  It can't always do that, and we will
# only know for sure later when we call predict with the calculated type.
# In this routine we call plotmo_predict with type=NULL to get all the
# response columns.  The dots are passed on to predict.
#
# If you don't need the response, set get.y=FALSE to reduce the amount of processing.

plotmo_meta <- function(object, type, nresponse, trace,
                        avoid.predict=FALSE, residtype=type,
                        msg.if.predictions.not.numeric=NULL, ...)
{
    type      <- plotmo_type(object, trace, "plotmo", type, ...)
    residtype <- plotmo_residtype(object, trace, "plotmo", residtype, ...)
    assignInMyNamespace("trace.call.global", trace) # trace call to resids, etc
    if(avoid.predict) {
        trace2(trace,
            "\n----Metadata: plotmo_resids(object, type=\"%s\", nresponse=NULL)\n",
            type)
        plotmo_resids <- plotmo_resids(object, type, residtype,
                                       nresponse=NULL, trace, ...)$resids
        if(is.null(plotmo_resids)) {
            if(trace >= 1)
                printf("residuals() was unsuccessful\n")
            avoid.predict <- FALSE # fall back to using predict
        } else {
            # trace2(trace,
            #     "got residuals using residuals(object, type=\"%s\", ...)\n", type)
            # use fitted rather than predict (TODO not right but ok for plotres)
            trace2(trace, "\n----Metadata: plotmo_fitted with nresponse=NULL\n")
            # nresponse=NULL so this returns multiple columns if a mult respe model
            plotmo_fitted <- plotmo_fitted(object, trace, nresponse=NULL, type, ...)
            yhat <- plotmo_fitted$fitted
            if(!inherits(object, "earth"))
                colnames(fitted) <- NULL # ensure get.resp.name.from.metadata doesn't use this
        }
    }
    if(!avoid.predict) {
        trace2(trace,
               "\n----Metadata: plotmo_predict with nresponse=NULL and newdata=NULL\n")
        # newdata=3 for efficiency
        plotmo_predict <- plotmo_predict(object, newdata=3, nresponse=NULL,
                            type, expected.levs=NULL, trace, inverse.func=NULL, ...)
        yhat <- plotmo_predict$yhat
        if(!is.null(msg.if.predictions.not.numeric)) {
            if(!is.null(plotmo_predict$resp.levs))
                stopf("%s when the predicted response is a factor",
                      msg.if.predictions.not.numeric)
            if(plotmo_predict$resp.class[1] == "character")
                stopf("%s when the predicted values are strings",
                      msg.if.predictions.not.numeric)
        }
        trace2(trace, "\n----Metadata: plotmo_fitted with nresponse=NULL\n")
        # nresponse=NULL so this returns multiple columns if a multiple response model
        plotmo_fitted <- plotmo_fitted(object, trace, nresponse=NULL, type, ...)
    }
    assignInMyNamespace("trace.call.global", 0)
    yfull <- NULL # plotmo_y with nresponse=NULL
    trace2(trace, "\n----Metadata: plotmo_y with nresponse=NULL\n")
    # nresponse=NULL so this returns multiple columns if a multi response model
    yfull <- plotmo_y(object, nresponse=NULL, trace,
                      expected.len=nrow(plotmo_fitted$fitted))$y
    nresponse.org <- nresponse

    nresponse <- plotmo_nresponse(yhat, object, nresponse, trace,
                    sprintf("predict.%s", class(object)[1]), type)
    stopifnot(!is.na(nresponse))
    trace2(trace,
        "nresponse=%g%s ncol(fitted) %d ncol(predict) %d ncol(y) %s\n",
        nresponse,
        if(identical(nresponse, nresponse.org))
            ""
        else
            sprintf(" (was %s)",
                if(is.character(nresponse.org)) paste0("\"", nresponse.org, "\"")
                else                            paste(nresponse.org)),
        NCOL(plotmo_fitted$fitted), NCOL(predict),
        sprintf("%d", NCOL(yfull)))

    y.as.numeric.mat <- NULL # y as single column numeric mat, only the nresponse column
    nresponse.y <- nresponse

    trace2(trace, "\n----Metadata: plotmo_y with nresponse=%g\n", nresponse)
    if(ncol(yfull) == 1 && nresponse.y > 1) {
        # e.g. lda(survived~., data=etitanic) with predict(..., type="post")
        nresponse.y <- 1
        trace2(trace,
            "nresponse=%d but for plotmo_y using nresponse=1 because ncol(y) == 1\n",
            nresponse)
    }
    assignInMyNamespace("trace.call.global", trace) # trace how we get the response
    y.as.numeric.mat <-
        plotmo_y(object, nresponse.y, trace, nrow(plotmo_fitted$fitted))$y
    assignInMyNamespace("trace.call.global", 0)

    resp.name <- get.resp.name.from.metadata(nresponse, trace,
                    yhat, plotmo_fitted$fitted, yfull, nresponse.y)

    resp.levs <- plotmo_resplevs(object, plotmo_fitted, yfull, trace)

    trace2(trace, "\n----Metadata: done\n\n")

    fitted <- plotmo_fitted$fitted

    list(
      yfull            = yfull,            # response as a data.frame, all columns

      y.as.numeric.mat = y.as.numeric.mat, # response as a single col numeric mat
                                           # only the nresponse column

      fitted           = fitted,    # fitted response as a data.frame (all columns)

      type             = type,      # type for predict()
                                    # always a string (converted from NULL if necesssary)

      residtype        = residtype, # type for residuals()
                                    # always a string (converted from NULL if necesssary)

      nresponse        = nresponse, # col index in the response (converted from NA if necessary)

      resp.name        = resp.name, # our best guess for the response name (may be NULL)

      resp.levs        = resp.levs) # levels of y before conversion to numeric (may be NULL)
                                    # necessary to convert predicted strings to factors
}
get.resp.name.from.metadata <- function(nresponse, trace,
                                        yhat, fitted, yfull, nresponse.y)
{
    # the order we look for the response name below seems to work but is not cast in stone
    if(is.factor(yhat[,1])) {
        # this prevents us putting a misleading first level name in plot headings
        resp.name <- NULL
        trace2(trace, "response name is NULL because is.factor(yhat[,1])\n")
    } else if(!is.null(colnames(yhat)) && nresponse <= length(colnames(yhat))) {
        resp.name <- colnames(yhat)[nresponse]
        trace2(trace, "got response name \"%s\" from yhat\n", resp.name)
    } else if(!is.null(yfull) && !is.null(colnames(yfull))) {
        resp.name <- colnames(yfull)[nresponse.y]
        trace2(trace, "got response name \"%s\" from yfull\n", resp.name)
    } else if(nresponse < length(colnames(fitted))) {
        resp.name <- colnames(fitted)[nresponse]
        trace2(trace, "got response name \"%s\" from plotmo_fitted\n", resp.name)
    } else {
        resp.name <- NULL
        trace2(trace, "response name is NULL\n")
    }
    resp.name
}
# Init resp.levs (the factor levels of the original response, may be NULL).
# The resp.levs is used if predict() returns strings (and therefore
# we must convert to them to a factor with the correct levels).

plotmo_resplevs <- function(object, plotmo_fitted, yfull, trace)
{
    levels.yfull <-
        if(is.null(yfull))
            NULL
        else if(length(dim(yfull)) == 2)
            levels(yfull[,1])
        else
            levels(yfull[1])

    if(!is.null(object[["levels"]])) {
        resp.levs <- object[["levels"]] # levels stored with earth
        trace2(trace, "got resp.levs from object$levels\n")
    } else if(!is.null(levels.yfull)) {
        resp.levs <- levels.yfull
        trace2(trace, "got resp.levs from yfull\n")
    } else if(!is.null(plotmo_fitted$resp.levs)) {
        resp.levs <- plotmo_fitted$resp.levs
        trace2(trace, "got resp.levs from plotmo_fitted$resp.levs\n")
    } else {
        resp.levs <- NULL
        trace2(trace, "resp.levs is NULL\n")
    }
    if(trace >= 2 && !is.null(resp.levs))
        printf("response levels: %s\n", paste.trunc(resp.levs))
    resp.levs
}
# This is used for processing "model response" variables such as the
# return value of predict(), fitted(), and residuals().
#
#
# If nresponse=NULL, return a data.frame but with y otherwise unchanged.
#
# Else return a numeric 1 x n matrix (regardless of the original class of y).
#   If nresponse is an integer, return only the specified column.
#   If nresponse=NA, try to convert it to a column index, error if cannot
#
# If subset is not NULL, take the specified subset of rows.
#
# If !is.null(nresponse) and y is character vector then convert it to a factor.
# expected.levs is used to do this (and not for anything else).
#
# returns list(y, resp.levs, resp.class)

process.y <- function(y, object, type, nresponse,
                      expected.len, expected.levs, trace, fname)
{
    if(is.null(y))
        stop0(fname, " NULL")
    if(length(y) == 0)
        stop0(fname, " zero length")
    print_summary(y, sprintf("%s returned", fname), trace)
    if(is.list(y) && !is.data.frame(y)) # data.frames are lists, hence must check both
        stop0(fname, " list, was expecting a vector, matrix, or data.frame\n",
              "       list(", list.as.char(y), ")")
    returned.resp.levs <- if(length(dim(y)) == 2) levels(y[,1]) else levels(y[1])
    resp.class <- class(y[1])
    colnames <- NULL
    resp.name <- NA
    dimy <- dim(y)
    if(length(dimy) == 3 && dimy[3] == 1) # hack for glmnet multnet objects
        y <- y[,,1]
    if(is.null(nresponse))
        y <- my.data.frame(y, trace, stringsAsFactors=FALSE)
    else {
        check.integer.scalar(nresponse, min=1, na.ok=TRUE, logical.ok=FALSE, char.ok=TRUE)
        nresponse <- plotmo_nresponse(y, object, nresponse, trace, fname, type)
        stopifnot(!is.na(nresponse), nresponse >= 1, nresponse <= NCOL(y))
        resp.name <- colname(y, nresponse, fname)
        y <- get.specified.col.and.force.numeric(y, nresponse, resp.name,
                                                 expected.levs, trace, fname)
        if(!is.na(nresponse) && nresponse > 1)
            print_summary(y, sprintf("%s returned", fname), trace,
                sprintf(" after selecting nresponse=%d", nresponse))
    }
    any.nas <- anyNA(y)
    any.non.finites <- FALSE
    # we use apply below because is.finite doesn't work for dataframes
    any.non.finites <- !any.nas &&
                       any(apply(y, 2, function(y) is.numeric(y) && !is.finite(y)))
    if(any.nas) {
        trace2(trace, "\n")
        warning0("NAs returned by ", fname)
    }
    if(any.non.finites) {
        trace2(trace, "\n")
        warning0("non-finite values returned by ", fname)
    }
    # Error message for the aftermath of:
    #   "Warning: 'newdata' had 100 rows but variable(s) found have 30 rows"
    if(!is.null(expected.len) && expected.len != nrow(y))
        stopf("%s returned the wrong length (got %d but expected %d)",
               fname[1], nrow(y), expected.len[1])

    print_summary(y,
                  sprintf("%s after processing with nresponse=%s is ",
                          fname,
                          if(is.null(nresponse)) "NULL" else format(nresponse)),
                  trace)
    list(y          = y, # n x 1 numeric, column name is original y column name
         resp.levs  = returned.resp.levs,
         resp.class = resp.class)
}
# always returns a one column numeric matrix
get.specified.col.and.force.numeric <- function(y, nresponse, resp.name,
                                                expected.levs, trace, fname)
{
    # nresponse=NA is not allowed at this point
    stopifnot(is.numeric(nresponse), length(nresponse) == 1, !is.na(nresponse))
    if(length(dim(y)) == 2)
        y <- y[, nresponse]
    else
        stopifnot(nresponse == 1)
    if(is.factor(y[1])) {
        trace2(trace, "converted to numeric from factor with levels %s\n",
               quotify.trunc(levels(y)))
        # TODO this may be a bogus warning
        if(!is.null(expected.levs) && !all.equal(expected.levs, levels(y[1])))
            warning0(fname, " returned a factor with levels ",
                     quotify.trunc(levels(y[1])),
                     " (expected levels ", quotify.trunc(expected.levs), ")")
    } else if(is.character(y[1])) {   # convert strings to factor
        old.y <- y
        y <- if(is.null(expected.levs)) factor(y)
             else                       factor(y, levels=expected.levs)
        trace2(trace, "converted to numeric from strings using factor levels %s\n",
               quotify.trunc(expected.levs))
        which <- (1:length(y))[is.na(y)]
        if(length(which)) {
            cat("\n")
            print_summary(old.y, fname, trace=2)
            cat("\n")
            printf("%s[%d] was %s and was converted to \"%s\"\n",
                fname, which[1], old.y[which[1]],
                if(is.na(y[which[1]])) "NA" else paste0("\"", y[which[1]], "\""))
            cat("\n")
            stopf("could not convert strings returned by %s to a factor (see above)",
                   fname)
        }
    }
    if(any(!is.double(y))) # convert logical or factor to double
        y <- as.vector(y, mode="numeric")
    y <- as.matrix(y)
    colnames(y) <- resp.name
    y
}
plotmo_nresponse <- function(y, object, nresponse, trace, fname, type="response")
{
    check.integer.scalar(nresponse, min=1, na.ok=TRUE, logical.ok=FALSE, char.ok=TRUE)
    colnames <- safe.colnames(y)
    nresponse.org <- nresponse
    if(is.na(nresponse)) {
        nresponse <- plotmo.convert.na.nresponse(object, nresponse, y, type)
        if(!is.na(nresponse)) {
            if(trace > 0 && nresponse != 1)
                printf("set nresponse=%s\n", paste(nresponse))
        } else { # nresponse is NA
            # fname returned multiple columns (see above) but nresponse is not specified
            cat("\n")
            print_summary(y, fname, trace=2)
            cat("\n")
            colnames <- NULL
            if(is.null(colnames) && !is.null(dim(y)))
                colnames <- colnames(y)
            icol <- min(3, NCOL(y))
            if(is.null(colnames))
                msg1 <- sprintf("%s\n       Example: nresponse=%d",
                    "Use the nresponse argument to specify a column.",
                    icol)
            else
                msg1 <- sprintf(
                    "%s\n          Example: nresponse=%s\n          Example: nresponse=%d",
                    "Use the nresponse argument to specify a column.",
                    quotify(if(is.na(colnames(y)[icol])) colname(y, 1) else colname(y, icol)),
                    icol)
            stopf(
"%s returned multiple columns (see above) but nresponse is not specified\n       %s",
                  fname, msg1)
        }
    } else if (is.character(nresponse)) {
        # convert column name to column index
        stopifnot.string(nresponse)
        if(is.vector(y))
            stop0("nresponse=\"", nresponse,
                  "\" cannot be used because the predicted response is a vector (it has no columns)")
        if(is.factor(y))
            stop0("nresponse=\"", nresponse,
                  "\" cannot be used because the predicted response is a factor (it has no columns)")
        if(is.null(colnames))
            stop0("nresponse=\"", nresponse,
                  "\" cannot be used because the predicted response has no column names")
        # TODO investigate [1] e.g. for plotmo(a1h.update2, nresponse="numd")
        nresponse <- imatch.choices(nresponse, colnames, err.msg.has.index=TRUE)[1]
    }
    check.integer.scalar(nresponse, min=1, na.ok=TRUE, logical.ok=FALSE, char.ok=TRUE)
    # note that msg is inhibited for trace<0, see trace1 in plotmo_rinfo
    if(nresponse > NCOL(y) && trace >= 0) {
        cat("\n")
        print_summary(y, fname, trace=2)
        cat("\n")
        check.index(nresponse, "nresponse", y,
                is.col.index=TRUE, allow.negatives=FALSE, treat.NA.as.one=TRUE)
    }
    if(trace >= 2 && (is.na(nresponse.org) || nresponse.org != nresponse))
        cat0("converted nresponse=",
             if(is.character(nresponse.org))
                paste0("\"", nresponse.org, "\"") else nresponse.org,
             " to nresponse=", nresponse, "\n")
    nresponse
}
plotmo.convert.na.nresponse <- function(object, nresponse, yhat, type="response")
{
    UseMethod("plotmo.convert.na.nresponse")
}
plotmo.convert.na.nresponse.default <- function(object, nresponse, yhat, type)
{
    stopifnot(is.na(nresponse))
    if(NCOL(yhat) == 1)
        1
    else if(NCOL(yhat) == 2 && substr(type, 1, 1) == "p")
        2     # probability (also works for posterior as in lda models)
    else
        NA
}
