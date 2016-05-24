# response.R: plotmo functions to get the response column from the given newdata
#                    mostly used for calculating RSq on newdata
#
# TODO overall structure here needs a bit of work

plotmo_rsq <- function(object, newdata=NULL, trace=0, nresponse=NA, type=NULL, ...)
{
    init.global.data() # needed if plotmo has never been invoked
    object.name <- quote.deparse(substitute(object))
    # TODO revisit, not really reliable because it may use parent.frame
    attr(object, ".Environment") <- get.model.env(object, object.name, trace)
    meta <- plotmo_meta(object, type, nresponse, trace, ...)
    plotmo_rsq1(object=object, newdata=newdata, trace=trace, meta=meta, ...)
}
plotmo_rsq1 <- function(object, newdata, trace, meta, ...)
{
    trace2(trace, "--plotmo_response for plotmo_rsq1\n")
    ynew <- plotmo_response(object=object, newdata=newdata, trace=max(0, trace),
                            nresponse=meta$nresponse, type=meta$type,
                            meta=meta, ...)
    trace2(trace, "--plotmo_predict for plotmo_rsq1\n")
    yhat <- plotmo_predict(object=object, newdata=newdata, nresponse=meta$nresponse,
                           type=meta$type, expected.levs=meta$expected.levs,
                           trace=trace, inverse.func=NULL, ...)$yhat
    if(ncol(yhat) != 1 || ncol(ynew) != 1 || nrow(yhat) != nrow(ynew)) {
        if(trace > -1) {
            printf("\n")
            print_summary(ynew, "response", trace=2)
            printf("\n")
            print_summary(yhat, "predicted values", trace=2)
            printf("\n")
        }
        stop0("response or predicted values have the wrong dimensions%s",
              if(trace > -1) " (see above)" else "")
    }
    get.weighted.rsq(ynew, yhat)
}
# If newdata is null, return the fitted response (same as plotmo_y).
#
# Else extract the response column from newdata.
# Use the model object to figure out which column is the response column.

plotmo_response <- function(object, newdata=NULL, trace=0,
                            nresponse=NA, type=NULL, meta=NULL, ...)
{
    print_summary(newdata, "--plotmo_response for newdata", trace)
    object.name <- quote.deparse(substitute(object))
    # TODO revisit, not really reliable because it may use parent.frame
    attr(object, ".Environment") <- get.model.env(object, object.name, trace)
    if(is.null(meta))
        meta <- plotmo_meta(object, type, nresponse, trace,
                            msg.if.predictions.not.numeric="RSq is not available",
                            ...)
    expected.len <- if(is.null(newdata)) NROW(meta$fitted) else NROW(newdata)
    if(is.null(newdata))
        y <- plotmo_y(object, meta$nresponse, trace,
                      expected.len=expected.len, resp.levs=meta$resp.levs)$y
    else if(length(dim(newdata)) != 2)
        stop0("plotmo_response: newdata must be a matrix or data.frame")
    else if(is.null(object$terms))
        y <- response.from.xy.model(object, newdata, trace, meta$resp.name)
    else # model has terms, presumably it was created with a formula
        y <- get.x.or.y.from.model.frame(object, field="y", trace, naked=FALSE,
                                         na.action=na.pass, newdata)$x
    if(!is.good.data(y, "response", trace, check.colnames=FALSE))
        stop0("response with newdata", format.err.field(y, "response", trace))
    y <- cleanup.x.or.y(object, y, "y", trace, check.naked=FALSE)
    if(!is.good.data(y, check.colnames=FALSE))
        stop0("response with newdata", format.err.field(y, "response", trace))
    y <- convert.glm.response(object, y, trace) # TODO test this and factor responses
    # TODO following will sometimes give the wrong results?
    if(!is.null(meta$nresponse) && meta$nresponse > NCOL(y)) {
        # if(trace >= 1)
        #     printf("forcing nresponse=%g to 1 because response has one column\n",
        #            nresponse)
        meta$nresponse <- 1
    }
    process.y(y, object, meta$type, meta$nresponse, expected.len=expected.len,
              meta$resp.levs, trace, "plotmo_response")$y
}
# the model was created with the x,y interface (no formula)

response.from.xy.model <- function(object, newdata, trace, resp.name)
{
    stopifnot(is.null(object$terms)) # shouldn't be here if model has formula
    if(!is.character(resp.name) || length(resp.name) != 1 || !nzchar(resp.name)) {
        if(trace > 2) {
            printf("\nresp.name:\n")
            print(resp.name)
            printf("\n")
        }
        stop0("could not get the response name")
    }
    trace2(trace, "response.from.xy.model: resp.name \"%s\"\n", resp.name)
    # following is for e.g. trees$Volume to Volume in earth(trees[,1:2], trees$Volume)
    resp.name <- sub(".*\\$", "", resp.name)
    # Hackery: look for responses of the form trees[,3] or trees[,3,drop=FALSE]
    # This happens if you build a model like lm(trees[,1:2], trees[,3])
    if(grepl("\\[.*,.+\\]", resp.name)) {
        col.name <- sub("[^,]*,", "", resp.name) # delete up to the comma and the comma
        col.name <- gsub(",.*",   "", col.name)  # delete (2nd) comma if any, and all after
        col.name <- gsub("\\]",   "", col.name)  # delete final ] if above gsub didn't do it
        # print a message because we don't always get this right
        if(trace >= 0)
            printf("Assuming response %s implies that the response column is %s\n",
                   resp.name, paste(col.name))
        # the following will do something like eval(3, env)
        col.index <- try.eval(parse(text=col.name), model.env(object), trace=trace, expr.name=col.name)
        if(is.try.err(col.index))
            stopf("could not parse the response name %s", resp.name)
        if(is.null(colnames(newdata)))
            resp.name <- paste0("newdata[,", col.index, "]")
        else # TODO is the following correct?
            resp.name <- paste0(colnames(newdata)[col.index])
        y <- newdata[, col.index, drop=FALSE]
    } else { # resp.name doesn't have [] in it, hopefully it's just a name
        colnames.newdata <- colnames(newdata)
        if(is.null(colnames.newdata))
            stop0("cannot get response from newdata because newdata has no column names")
        which <- which(colnames.newdata == resp.name)
        if(length(which) == 0)
            stop0("no column names in newdata match the original response name\n",
                  sprintf("       Response name: %s\n", resp.name),
                  "       Column names in newdata: ", paste.collapse(colnames.newdata))
        if(length(which) > 1)
            stopf("multiple column names in newdata match the original response name %s",
                  resp.name)
        y <- newdata[, colnames.newdata[which], drop=FALSE]
    }
    y
}
