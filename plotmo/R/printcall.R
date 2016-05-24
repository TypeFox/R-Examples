# printcall.R: functions for printing call information

# If call is specified, print it (where call is from match.call or similar).
# Else use the call stack to determine the call. The n arg tells us how
# far to go back in the call stack.
#
# Examples: printcall()     describe the call to the current function
#           printcall(n=2)  describe the call to the caller of the current function
#           printcall(call) describe call where call is from match.call or similar

printcall <- function(prefix="", call=NULL, all=FALSE, n=1)
{
    # check prefix and n here, other args checked in call.as.char
    stopifnot.string(prefix, allow.empty=TRUE)
    stopifnot(is.numeric(n))
    call <- call.as.char(call, all, n+1)
    printf.wrap("%s%s\n", prefix, call)
}
# returns args and concise description of their values, dots are included
# all=TRUE to include all formal args (not always avail e.g. for primitives)
#
# TODO Does not expand the dots (just prints "..."), need fixed version of match.call
#      To expand the dots see e.g. higher.call.to.deprefix (but that would only work
#      here if dots for caller at n where the same as the dots to printcall

call.as.char <- function(call=NULL, all=FALSE, n=1)
{
    stopifnot(is.numeric(all) || is.logical(all), length(all) == 1)
    stopifnot(is.numeric(n), length(n) == 1, n > 0)
    if(is.null(call))
        call <- match.call2(all=all, n=n+1) # +1 to skip call to call.as.char
    else if(all) # we have the call but not the func itself, so can't get formals
        stop("all=TRUE is not allowed when the call argument is used")
    fname <- fname.from.call(call)
    if(all) {
        formals <- formals(attr(call, "sys.function"))
        call[[1]] <- NULL               # delete func name from call, leave args
        formals[["..."]] <- NULL        # delete ... in formal args if any
        call <- merge.list(formals, call)
    } else
        call[[1]] <- NULL               # delete func name from call, leave args
    rv <- paste(fname, "(", list.as.char(call, maxlen=50), ")", sep="")
    attr(rv, "fname") <- fname # needed for alignment with nchar in printcall
    rv
}
# Similar to match.call but with args "all" and "n".
# Also, this always returns a call, even if it is merely "unknown()".
# So you can safely call it with any n (although n must be a positive int).

match.call2 <- function(all=FALSE, n=1)
{
    stopifnot(is.numeric(all) || is.logical(all), length(all) == 1)
    stopifnot(is.numeric(n), length(n) == 1, n > 0)
    # get sys.function and sys.call for the given n, needed for match.call
    sys.function <- try(sys.function(-n), silent=TRUE)
    if(is.try.err(sys.function) || is.null(sys.function)) # typically "not that many frames"
        return(call("unknown"))
    sys.call <- try(sys.call(-n), silent=TRUE)
    if(is.try.err(sys.call) || is.null(sys.call))
        return(call("unknown"))
    # TODO following can cause incorrect "... used in a situation where it does not exist"
    #      R version 3.1.4 will fix that issue in match.call (I hope)
    # envir <- parent.frame(n+1) # use when new version of match.call is ready
    call <- try(match.call(definition=sys.function, call=sys.call, expand.dots=TRUE),
                silent=TRUE)
    if(is.try.err(call)) {
        # match.call failed, fallback to a weaker description of call
        # no expansion of dots and no arg values :(
        call <- sys.call
    }
    attr(call, "sys.function") <- sys.function
    call
}
callers.name <- function(n=1)
{
    stopifnot(is.numeric(n), length(n) == 1, floor(n) == n, n >= 0)
    call <- try(sys.call(-(n+1)), silent=TRUE)
    fname.from.call(call) # will also check if try error
}
fname.from.call <- function(call) # call was obtained using sys.call() or similar
{
    if(is.try.err(call))
        return("unknown") # most likely n was misspecified (too big)
    if(is.null(call)) # e.g. NULL->source->withVisible->eval->eval->print->test->callers.name
        return("NULL")
    caller <- as.list(call)[[1]]
    if(is.name(caller))   # e.g. foo3(x=1)
        caller <- as.character(caller)
    else {                # class(caller) is "call" e.g. plotmo::localfunc(x=1)
        stopifnot(is.call(call))
        caller <- format(caller)
    }
    if(grepl("function (", substr(caller[1], 1, 10), fixed=TRUE))
        paste0("function(", paste.trunc(strip.space(substring(caller, 11))), ")")
    else
        paste.trunc(strip.space(caller))
}
# if EVAL is FALSE this will print something like xlim=..1, ylim=..2
# TODO add n arg when match.call is fixed (R version 3.2.0)
# TODO also then make this callable as printdots() instead of printdots(...)

printdots <- function(..., EVAL=TRUE, PREFIX=sprintf("%s dots: ", callers.name))
{
    sys.call <- as.list(sys.call())
    ensure.dots.present(sys.call)
    callers.name <- callers.name()
    printf.wrap("%s%s\n", PREFIX, dots.as.char(..., EVAL=EVAL))
}
dots.as.char <- function(..., EVAL=TRUE)
{
    sys.call <- as.list(sys.call())
    ensure.dots.present(sys.call)
    dots <- match.call(expand.dots=FALSE)$...
    if(is.null(dots))
        return("no dots")
    if(EVAL) {
        stopifnot(is.numeric(EVAL) || is.logical(EVAL), length(EVAL) == 1)
        dots <- eval.dotlist(dots)
    }
    list.as.char(dots)
}
# issue error message if ... wasn't used in the call to dots.as.char
ensure.dots.present <- function(sys.call)
{
    dots.present <- FALSE
    for(i in seq_len(length(sys.call)))
        if(sys.call[i] == "...")
            dots.present <- TRUE
    if(!dots.present)
        stop0("dots.as.char should be invoked with dots, for example dots.as.char(...)")
}
