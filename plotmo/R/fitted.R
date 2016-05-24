# fitted.R: plotmo functions for getting the fitted data for an arbitrary model

# Like fitted() but will get fitted response even if not already with object.
# Returns an n x 1 matrix (unless nresponse=NULL then returns an n x q dataframe?).
# The returned columns may not be named.
# The type and dots args are used if the call to fitted(object) fails.

plotmo_fitted <- function(object, trace, nresponse, type, ...)
{
    if(!is.null(nresponse))
        check.integer.scalar(nresponse, min=1, na.ok=TRUE, logical.ok=FALSE, char.ok=TRUE)
    fitted <- try(call.dots(stats::fitted, DROP="*", KEEP="PREFIX",
                    # following prevents reprint of fitted msg if fail
                    TRACE=if(trace <= 0) -1 else (trace >= 2 || trace.call.global),
                    force.object=object, ...),
                  silent=trace <= 1)
    if(!is.try.err(fitted) && !is.null(fitted)) {
        # if(trace.call.global >= 1 && trace < 2)
        #     print_summary(fitted, "fitted is ", details=-1)
        temp <- process.y(fitted, object, type, nresponse, expected.len=NULL,
                          expected.levs=NROW(fitted), trace, "fitted(object)")
        fitted <- temp$y
    } else { # fitted(object) failed
        if(trace >= 1)
            printf("fitted() was unsuccessful, will use predict()\n")
        type <- plotmo_type(object, trace, "plotmo", type, ...)
        # we have already printed call to predict so clear trace flag
        # (this is dependent on the sequence of calls in plotmo_meta)
        assignInMyNamespace("trace.call.global", 0)
        temp <- plotmo_predict(object, newdata=NULL, nresponse, type,
                    expected.levs=NULL, trace=trace, inverse.func=NULL, ...)
        fitted <- temp$yhat
        trace2(trace, "got fitted values by calling predict (see above)\n")
    }
    if(!is.null(colnames(fitted)))
        colnames(fitted) <- sub(".*\\$", "", colnames(fitted)) # trees$Volume to Volume
    list(fitted    = fitted, # n x 1 numeric unless nresponse=NULL
         resp.levs = temp$resp.levs)
}
