# prolog.R: plotmo.prolog functions, called at the start of plotmo and plotres

plotmo.prolog <- function(object, object.name, trace, ...) # gets called at the start of plotmo
{
    UseMethod("plotmo.prolog")
}
plotmo.prolog.default <- function(object, object.name, trace, ...)
{
    # prevent confusing downstream errors by doing an initial check here
    if(is.null(object$call) && is.null(object[["x"]]))
        stopf("%s does not have a 'call' field or %s",
              object.name,
              if(is.null(object[["y"]])) "'x' and 'y' fields"
              else                       "an 'x' field")
    object
}
