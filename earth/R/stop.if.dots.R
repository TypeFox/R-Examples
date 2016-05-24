# stop.if.dots.R:

stop.if.dots <- function(...)
{
    dots <- match.call(expand.dots=FALSE)$...
    if(length(dots))
        dots.used.err(STOPFUNC=base::stop, MSG=": unrecognized", ...)
}

warn.if.dots <- function(...)
{
    dots <- match.call(expand.dots=FALSE)$...
    if(length(dots))
        dots.used.err(STOPFUNC=base::warning, MSG=" ignored", ...)
}
dots.used.err <- function(..., STOPFUNC, MSG)
{
    callers.name <- callers.name(n=2)
    dots <- match.call(expand.dots=FALSE)$...
    for(idot in seq_along(dots)) # STOPFUNC is either stop() or warning()
        STOPFUNC(callers.name, MSG, describe.dot(dots, idot), call.=FALSE)
}
describe.dot <- function(dots, idot)
{
    nchar <- nchar(names(dots)[idot])
    if(length(nchar) && nchar > 0)
        return(sprintf(" argument '%s'", names(dots[idot])))
    # the argument that was passed in dots is unnamed
    call <- call.as.char(n=5) # n=5 to describe call to caller of stop.if.dots
    sprintf(" unnamed argument\n       The call was %s",
            paste0(strwrap(call, width=max(40, getOption("width")-20), exdent=25),
                   collapse="\n"))
}
