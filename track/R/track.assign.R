track.assign <- function(x, value, pos=1, envir=as.environment(pos), flush=TRUE) {
    if (!is.character(x) || length(x)!=1)
        stop("x must be a character vector of length 1")
    if (!env.is.tracked(envir))
        stop("envir ", environmentName(envir), " is not tracked")
    if (length(tracked(list=x, envir=envir))==0)
        track(list=x, envir=envir)
    assign(x, value, envir=envir)
    if (flush)
        track.flush(list=x, envir=envir)
    invisible(NULL)
}
