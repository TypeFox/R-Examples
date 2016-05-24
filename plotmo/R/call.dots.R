# call.dots.R: functions to handle prefixed dot arguments

call.dots <- function(
    FUNC     = NULL,
    ...,
    PREFIX   = NULL,
    DROP     = "*",
    KEEP     = "PREFIX",
    TRACE    = 0,
    FNAME    = if(is.character(FUNC)) FUNC
               else trunc.deparse(substitute(FUNC)),
    FORMALS  = NULL,  # formal args of FUNC (needed because CRAN doesn't allow :::)
    SCALAR   = FALSE, # see argument "scalar" in eval.dotlist
    CALLARGS = NULL,
    CALLER   = NULL)
{
    stopifnot(is.logical(TRACE) || is.numeric(TRACE), length(TRACE) == 1)
    TRACE <- as.numeric(TRACE)
    if(TRACE >= 2) {
        if(is.null(CALLER))
            CALLER <- callers.name()
        printf("%s invoked call.dots\n", CALLER)
    }
    if(is.null(CALLARGS))
        CALLARGS <- callargs(call.dots)

    args <- deprefix(FUNC=FUNC, PREFIX=PREFIX, ..., DROP=DROP, KEEP=KEEP,
                     TRACE=TRACE, FNAME=FNAME, FORMALS=FORMALS,
                     SCALAR=SCALAR, CALLARGS=CALLARGS)

    do.call.trace(FUNC, args, FNAME, trace=TRACE)
}
call.plot <- function(
    FUNC    = NULL,
    ...,
    PREFIX  = NULL,
    TRACE   = 0,
    FORMALS = NULL,
    SCALAR  = FALSE)
{
    fname <- trunc.deparse(substitute(FUNC))
    callargs <- callargs(call.plot)
    caller <- callers.name() # function that invoked call.plot

    call.dots(FUNC=FUNC, PREFIX=PREFIX, ...,
              DROP="*",                # drop everything
              KEEP="PREFIX,PLOT.ARGS", # except args matching PREFIX and PLOT.ARGS
              TRACE=TRACE, FNAME=fname,
              FORMALS=FORMALS, SCALAR=SCALAR, CALLARGS=callargs, CALLER=caller)
}
deprefix <- function(
    FUNC     = NULL,
    ...,
    PREFIX   = NULL,
    DROP     = NULL,
    KEEP     = NULL,
    TRACE    = 0,
    FNAME    = if(is.character(FUNC)) FUNC
               else trunc.deparse(substitute(FUNC)),
    FORMALS  = NULL,
    SCALAR   = FALSE,
    CALLARGS = NULL)
{
    stopifnot(is.logical(TRACE) || is.numeric(TRACE), length(TRACE) == 1)
    TRACE <- as.numeric(TRACE)
    if(!is.null(FUNC))
        match.fun(FUNC) # check that FUNC is available and is a function
    FNAME <- init.fname(FNAME, FUNC, TRACE)
    higher.caller <- higher.caller.to.deprefix(..., FNAME=FNAME)
    PREFIX <- init.prefix(PREFIX, FUNC, FNAME)
    if(is.null(CALLARGS))
        CALLARGS <- callargs(deprefix)
    DROP <- expand.drop(DROP, PREFIX, FUNC, FORMALS)
    KEEP <- expand.drop(KEEP, PREFIX, FUNC, FORMALS, namedrop="KEEP",
                        callargs=CALLARGS, include.standard.prefixes=TRUE)
    dots <- match.call(expand.dots=FALSE)$...
    trace.prolog(TRACE, PREFIX, DROP, KEEP, dots, higher.caller)
    stopif.unamed.dot(dots, higher.caller, ...)
    org.dots <- dots
    if(!is.null(DROP))
        dots[grep(DROP, names(dots))] <- NULL
    stopifnot(!is.null(KEEP))
    for(name in names(org.dots))
        if(grepl(KEEP, name))
            dots[[name]] <- org.dots[[name]]
    trace.after.dropkeep(TRACE, dots)
    args <- deprefix.aux(FUNC, dots, PREFIX, FNAME, FORMALS, TRACE) # workhorse
    eval.dotlist(args, n=2, scalar=SCALAR) # n=2 for caller of deprefix e.g. call.dots
}
deprefix.aux <- function(func, dots, prefix, fname, formals, trace) # workhorse
{
    force <- "^force\\."    # "force." as a regex
    def   <- "^def\\."      # "def." as a regex

    # change prefix to a regex, "plot." becomes "^plot\."
    prefix <- paste0("^", gsub(".", "\\.", prefix, fixed=TRUE))

    groups <- list() # list with three elements: force, prefix, def args
    for(pref in c(force, prefix, def)) {
        # put args matching pref into group, with the prefix pre removed
        which <- grep(pref, names(dots))        # select only args matching pref
        group <- dots[which]                    # put them into the group
        group <- expand.dotnames(group, pref, func, fname, formals)
        names(group) <- sub(pref, "", names(group)) # remove prefix
        groups[[pref]] <- group
        dots[which] <- NULL                     # remove args in this group from dots
    }
    # dots is now just those arguments which did not have a special prefix
    dots <- expand.dotnames(dots, prefix="", func, fname) # "" matches anything
    args <- groups[[def]]                       # "def." args lowest precedence
    args <- merge.list(args, dots)              # next come remaining dots
    args <- merge.list(args, groups[[prefix]])
    args <- merge.list(args, groups[[force]])   # "force." args overrule all others
    args <- drop.args.prefixed.with.drop(args)
    order.args(args, trace)
}
# Argument names for plot functions.  We exclude "overall" par() args like
# mfrow that shouldn't be included when calling functions like plot(),
# lines(), or text().
#
# If specified in a DROP or KEEP string, the actual argument must exactly
# match the PLOT.ARGS argument to be dropped or kept --- abreviated actual
# args won't be matched (otherwise we would match too much, e.g. an actual
# arg "s" would match "srt").

PLOT.ARGS <- c("add", "adj", "bty", "cex", "cex.axis", "cex.lab", "cex.main",
    "cex.sub", "col", "col.axis", "col.lab", "col.main", "col.sub",
    "crt", "family", "font", "font", "font.axis", "font.lab", "font.main",
    "font.sub", "lend", "ljoin", "lmitre", "lty", "lwd", "main", "pch",
    "srt", "xaxp", "xaxs", "xaxt", "xlab", "xlim", "xlog", "xpd", "yaxp",
    "yaxs", "yaxt", "ylab", "ylim", "ylog")

# Arguments for par().  This list includes all par arguments except
# readonly arguments (e.g. cin) and unimplemented arguments (e.g. err).
# The actual argname must be an exact match to be recognized (no abbreviations).
# Following omitted because they change too much: col, lwd

PAR.ARGS <- c("adj", "ann", "ask", "bg", "bty", "cex", "cex.axis", "cex.lab",
    "cex.main", "cex.sub", "col.axis", "col.lab", "col.main",
    "col.sub", "crt", "err", "family", "fg", "fig", "fin", "font",
    "font.axis", "font.lab", "font.main", "font.sub", "lab", "las", "lend",
    "lheight", "ljoin", "lmitre", "lty", "mai", "mar", "mex",
    "mfcol", "mfg", "mfrow", "mgp", "mkh", "new", "oma", "omd", "omi",
    "pch", "pin", "plt", "ps", "pty", "srt", "tck", "tcl", "usr", "xaxp",
    "xaxs", "xaxt", "xlog", "xpd", "yaxp", "yaxs", "yaxt", "ylbias", "ylog")

# Arguments for par() which take a vector value (i.e. length of value is not one).

PAR.VEC <- c("fig", "fin", "lab", "mai", "mar", "mfcol", "mfg", "mfrow", "mgp",
             "oma", "omd", "omi", "pin", "plt", "usr", "xaxp", "yaxp")

# Arguments that are used for subplots in plotmo and similar programs.
#
# Useful for dropping all args that could conceivably be plotting
# arguments and will never(?) be a predict() or residuals() argument.
#
# When "PLOTMO.ARGS" is used in a DROP string, any actual arg _prefixed_
# with any of these is dropped (as opposed to PLOT.ARGS and PAR.ARGS we drop
# actual argnames that _exactly_ match argnames in PLOT.ARGS and PAR.ARGS).
#
# "nresiduals", is for back compat with old versions of plot.earth

PLOTMO.ARGS <- c( "caption.", "cex.", "col.", "contour.", "cum.",
    "degree1.", "degree2.", "density.", "filled.contour.", "font.", "func.",
    "grid.", "heatmap.", "image.", "jitter.", "legend.", "label.", "level.",
    "line.", "lines.", "lty.", "lty.", "lwd.", "main.", "mtext.",
    "nresiduals", "par.", "pch.", "persp.", "plot.", "plotmath.", "qq.",
    "qqline.", "pt.", "response.", "rug.", "smooth.", "text.", "title.",
    "vfont.")

# from now on in this module function defs are in alphabetic order

add.formals.to.drop <- function(drop, func, formals, namedrop)
{
    stopifnot(grepl("FORMALS", drop))
    if(is.null(func))
        stop0("\"FORMALS\" specified in ", namedrop, ", but FUNC is NULL")

    formals <- merge.formals(func, formals, must.exist=TRUE)
    formals <- paste0(formals, collapse=",") # vector to string

    drop <- sub("FORMALS[,]", "", drop)     # remove "FORMALS," from drop
    paste.drop(">FORMALS", formals, drop)   # add the formal args
}
# Return the names of the actual args passed to the caller of this function,
# ignoring args matching formals of the caller and ignoring dots.
#
# For example, for call.dots(foo, PREFIX="anything", x=1, y=1, ...), this
# function returns c("x", "y"), because x and y are in the argument list
# in the call to call.dots but don't match any of the formals of call.dots
# (as PREFIX does).  The "..." is ignored.
# TODO if these were forced we wouldn't need the force.argument
callargs <- function(func)
{
    # names of arguments passed to the func that invoked callargs
    # args passed in dots will not appear in names
    names <- names(sys.call(-1))
    names <- names[names != ""] # drop unnamed args

    # drop formal arguments (typically PREFIX, KEEP, etc.)
    names[!(names %in% names(formals(func)))]
}
# return string "a,b,c,d,e" if given c("a", "b,c", "d e")
# i.e. white space converted to comma, c() collapsed to single string
canonical.drop <- function(drop, namedrop)
{
    drop <- gsub(" +|,+", ",", drop)          # convert space or multi commas to comma
    drop <- gsub("^,+|,+$", "", drop)         # drop leading and trailing commas
    drop <- unlist(strsplit(drop, split=",")) # convert to a vector
    drop <- paste0(drop, collapse=",")        # collapse
    stopifnot.identifier(drop, namedrop, allow.specials=TRUE)
    drop
}
# TODO add this check elsewhere in earth and plotmo too
check.regex <- function(s) # check for some common regex errors
{
    if(grepl("||", s, fixed=TRUE))
        stop0("\"||\" in following regex matches everything:\n",
              "\"", s, "\"")
    if(grepl("^\\|", s))
        stop0("\"|\" at the start of the following regex matches everything:\n",
              "\"", s, "\"")
    if(grepl("\\|$", s))
        stop0("\"|\" at the end of the following regex matches everything:\n",
              "\"", s, "\"")
}
# convert drop to a regex, "x,y*,prefix." becomes "^x|^y.*|^prefix\."
convert.drop.to.regex <- function(drop)
{
    drop <- gsub(",", "|",   drop)             # change comma to |
    drop <- gsub(".", "\\.", drop, fixed=TRUE) # escape period, "plot." becomes "plot\."
    drop <- gsub("*", ".*",  drop, fixed=TRUE) # change * to .*

    # clean up, for example we now may have "||" in drop which must be changed to "|"
    for(iter in 1:2) { # two iterations seems sufficient in practice
        drop <- gsub(" +", "",  drop)             # delete spaces
        drop <- sub("^\\|", "", drop)             # delete | at at start
        drop <- sub("^\\|", "", drop)             # delete | at at end
        drop <- gsub("^^", "^", drop, fixed=TRUE) # change ^^ to single ^
        drop <- gsub("||", "|", drop, fixed=TRUE) # change || to |
    }
    # prepend ^ to match prefixes only, "x|y" becomes "^x|^y"
    drop <- unlist(strsplit(drop, split="|", fixed=TRUE))
    drop <- ifelse(substr(drop, 1, 1) == ">", drop, paste0("^", drop))
    drop <- paste0(drop, collapse="|")

    check.regex(drop) # sanity check for some common regex errors

    drop
}
# TODO add to test suite (although this is tested implicitly in the plotmo tests)
#      what happens if the argname is abbreviated and no formals to match against?
drop.args.prefixed.with.drop <- function(args)
{
    for(name in names(args)) if(grepl("^drop\\.", name)) {
        check.integer.scalar(args[[name]], logical.ok=FALSE, object.name=name)
        if(args[[name]] != 1)
            stop0(name, "=1 is not TRUE")
        args[[name]] <- NULL            # drop the drop.xxx argument itself
        name <- sub("drop.", "", name, fixed=TRUE) # delete "drop." from name
        # TODO allow dropping if just the prefix of name matches
        name <- paste0("^", name, "$")  # turn it into a regex for exact matching
        args[grep(name, names(args))] <- NULL  # drop args that exactly match name
    }
    args
}
# Only dot names that have the given prefix are considered.  Expand the
# suffix of each of those dot names to its full formal name using the
# standard R argument matching rules.
#
# Example: with prefix = "persp." and func = persp.default,
# "persp.sh" in dots gets expanded to "persp.shade", because
# "shade" is the full name of an argument of persp.default.
#
# Among other things, This makes it possible for deprefix to properly
# process two actual argument names that are different but both match
# the same formal argument name.
#
# It also helps prevent downstream name aliasing issues, because here we
# can pre-emptively check for argname matching problems, and issue clearer
# error messages than the standard R arg matching error messages.

expand.dotnames <- function(
    dots,
    prefix,         # a regex, not a plain string
    func = NULL,    # if NULL then we just check for duplicate args and go home
    fname,          # used only in error messages
    formals = NULL) # manual additions to the formal arg list of func
{
    stopifnot(is.list(dots))
    dot.names <- names(dots)
    matches <- grep(prefix, dot.names) # indices of arg which match prefix
    if(length(matches) == 0)
        return(list())
    if(is.null(func)) {
        duplicated <- which(duplicated(dot.names))
        if(length(duplicated))
            stop0("argument '", dot.names[duplicated[1]], "' for ",
                  fname, "() is duplicated")
        return(dots[matches])
    }
    # match against the formal arguments of func
    stopifnot(!is.null(dot.names))
    unexpanded.names <- dot.names
    formals <- merge.formals(func, formals)
    for(idot in matches) { # for all arguments which match prefix
        dot.name <- dot.names[idot]
        stopifnot(nzchar(dot.name))
        raw.prefix <- ""
        raw.dotname <- dot.name
        if(nzchar(prefix)) {
            # strip off the prefix substring in dot.name (we will put it back later)
            start <- regexpr(prefix, dot.name)
            stopifnot(start == 1) # prefix matches only prefixes
            stop <- start + attr(start, "match.length")
            stopifnot(stop > start)
            raw.prefix <- substr(dot.name, start=start, stop=stop-1) # as string not regex
            raw.dotname <- substring(dot.name, first=stop) # dotname with prefix removed
        }
        match <- charmatch(raw.dotname, formals)

        if(is.na(match)) {
            # No match, not necessarily a problem assuming FUNC has a dots formal arg.
            # We will allow FUNC to check for itself later (if someone calls it).
            NULL

        } else if(match == 0) { # multiple matches
            matches <- grep(paste0("^", raw.dotname), formals)
            stopifnot(length(matches) >= 2)
            stop0("'", raw.dotname, "' matches both the '", formals[matches[1]],
                  "' and '", formals[matches[2]], "' arguments of ", fname, "()")

        } else # single match, this is the ideal situation
            dot.names[idot] <- paste0(raw.prefix, formals[match]) # prepend prefix
    }
    stopifnot.expanded.dotnames.unique(dot.names, unexpanded.names,
                                       fname, formals, prefix)
    names(dots) <- dot.names
    dots
}
# returned the expanded the drop argument as a regex

expand.drop <- function(drop, prefix, func,
                formals=NULL, # manual additions to the formal arg list of func
                namedrop="DROP", callargs=NULL, include.standard.prefixes=FALSE)

{
    if(is.null(drop)) {
        if(include.standard.prefixes)
            return(paste0("^force.|^def.|^", prefix))
        else
            return(NULL)
    }
    drop <- canonical.drop(drop, namedrop)

    if(drop == "*")
        return(".*") # regex to match everything

    # TODO following is helpful in the trace print only if
    # you put special identifiers AFTER the other identifiers
    drop <- paste.drop(">EXPLICIT", drop, "")

    if(length(callargs) > 0)
        drop <- paste.drop(">CALLARGS,", paste0(callargs, "$", collapse=","), drop)

    if(include.standard.prefixes) {
        drop <- sub("PREFIX", "", drop) # delete "PREFIX" from drop, if present
        drop <- paste.drop(">PREFIX,", prefix, drop)
        drop <- paste.drop(">STANDARDPREFIXES,", "force.,def.,drop.", drop)
    } else
        drop <- paste.drop(">PREFIX,", sub("PREFIX", prefix, drop), "")

    if(grepl("FORMALS", drop))
        drop <- add.formals.to.drop(drop, func, formals, namedrop)

    temp <- paste.drop(">PLOT_ARGS,", paste0(PLOT.ARGS, "$", collapse=","), "")
    drop <- sub("PLOT.ARGS", temp, drop)

    temp <- paste.drop(">PAR_ARGS,", paste0(PAR.ARGS, "$", collapse=","), "")
    drop <- sub("PAR.ARGS", temp, drop)

    temp <- paste.drop(">PLOTMO_ARGS,", paste0(PLOTMO.ARGS, collapse=","), "")
    drop <- sub("PLOTMO.ARGS", temp, drop)

    convert.drop.to.regex(drop) # convert drop to a regex
}
higher.call.args <- function(..., CALLX, FNAME)
{
    stopifnot(is.list(CALLX))
    CALLX[1] <- NULL                  # remove fname from CALLX
    if(CALLX[length(CALLX)] == "...") # remove dots from CALLX
        CALLX[length(CALLX)] <- NULL
    args <- eval.dotlist(as.list(CALLX))
    # add dots to args, if they are not already in args
    dots <- as.list(match.call(expand.dots=FALSE)$...)
    arg.names <- names(args)
    dot.names <- names(dots)
    for(i in seq_along(dots)) {
        if(!(dot.names[i] %in% arg.names)) {
            list <- list(eval(dots[[i]]))
            names(list) <- dot.names[i]
            args <- append(args, list)
        }
    }
    args[[1]] <- as.name(FNAME)
    list.as.char(args)
}
# used only for tracing and error messages
# TODO simplify this and friends when match.call is working (R 3.2.0)
higher.caller.to.deprefix <- function(..., FNAME=FNAME)
{
    # search the stack looking for org caller of prefix e.g. call.plot
    sys.calls <- sys.calls()
    ncalls <- length(sys.calls)
    stopifnot(ncalls > 2)
    higher.fname <- "FUNC"
    try.was.used <- FALSE
    for(i in max(ncalls-10, 1) : ncalls) {
        fname <- paste(sys.calls[[i]])[1]
        if(grepl("^call\\.|^deprefix", fname))
            break
        if(grepl("^doTry|^try", fname))
            try.was.used <- TRUE
        else
            higher.fname <- fname
    }
    call  <- as.list(sys.calls[[i]])
    fname <- paste(call[[1]])
    if(try.was.used)
        higher.fname <- paste0(higher.fname, " via try ")
    # use try here for paranoia
    args <- try(higher.call.args(..., CALLX=call, FNAME=FNAME), silent=TRUE)
    if(is.try.err(args))
        args <- sprintf("%s, ...", FNAME)
    sprintf("%s called %s(%s)", higher.fname, fname, args)
}
init.fname <- function(FNAME, FUNC, TRACE)
{
    # check deparse(substitute(FUNC)) issued a good function name
    # e.g. FNAME will be "NULL" if FUNC is NULL
    if(is.null(FNAME) || length(FNAME) != 1 || FNAME == "NULL")
        FNAME <- "FUNC"
    stopifnot.string(FNAME)
    FNAME <- sub(".*:+", "", FNAME) # graphics::lines becomes lines
    stopifnot.identifier(FNAME, "FNAME")
    FNAME
}
init.prefix <- function(PREFIX, FUNC, FNAME)
{
    if(is.null(PREFIX)) {
        # automatic prefix, so check that we can generate it safely
        if(is.null(FUNC))
            stop0("PREFIX must be specified when FUNC is NULL")
        PREFIX <- sub("\\..*$", "", FNAME) # lines.default becomes lines
        # Was deprefix invoked using FUNC=FUNC or in a try block?
        # This won't catch all cases of FUNC=unusable.name but it helps
        # The stopifnot.identifier() below also helps.
        if(PREFIX %in% c("FUNC", "doTryCatch"))
            stop0("PREFIX must be specified in this context ",
                  "(because FNAME is \", fname, \")")
        PREFIX <- paste0(PREFIX, ".")      # add a period
        stopifnot.identifier(PREFIX, "the automatically generated PREFIX")
    }
    stopifnot.identifier(PREFIX, "PREFIX", allow.empty=TRUE)
    if(PREFIX == "")
        PREFIX <- ">NOPREFIX" # no argname can match this
    PREFIX
}
# return a char vector: formal() of func plus names in manform
# manform is manually specified formals
merge.formals <- function(func, manform, must.exist=FALSE)
{
    formals <- names(formals(func))
    if(!is.null(manform))
        formals <- c(formals,
                     strsplit(canonical.drop(manform, "manform"), ",")[[1]])
    if(must.exist) {
        if(length(formals) == 0)
            stop0("\"FORMALS\" specified but formals(FUNC) ",
                  "returned no formal arguments")
        if(length(formals[formals != "..."]) == 0)
            stop0("\"FORMALS\" specified but formals(FUNC) returned only \"...\"")
    }
    formals <- formals[formals != "..."]  # drop arg named "..." in formals, if any
    sapply(formals, stopifnot.identifier) # check that all names are valid
    unique(formals)
}
# Put the "anon" args first in the argument list.
# Then put args named "object", "x", etc. at the front of the list
# (after the anon args if any).  This is necessary because all the
# manipulation we have done has sadly done some reordering of the args
# (meaning that the order of the args supplied to call.dots is only
# partially retained).  The names object, x, etc. are usually what we want
# at the start for the predict and plot functions used with call.dots.

order.args <- function(args, trace)
{
    trace2(trace, "return dotnames      ")
    if(length(args)) {
        # order anonymous args on their names, then delete their names
        which <- which(grepl("^anon", names(args)))
        anon <- args[which]              # select args with "anon." prefix
        args[which] <- NULL              # remove them from the arg list
        anon <- anon[order(names(anon))] # order them on their names
        trace2(trace, "%s", paste0(names(anon), collapse=" "))
        names(anon) <- NULL              # delete their names
        args1 <- anon                    # anon args go first in the arg list
        # Put arguments named "object", "x", etc. first (after anon args if any).
        # We want mfrow and mfcol early so subsequent args like cex have the last say.
        for(argname in c("object", "x", "y", "type", "main",
                         "xlab", "ylab", "mfrow", "mfcol")) {
            args1[[argname]] <- args[[argname]]
            args[[argname]]  <- NULL # remove from args
        }
        args <- append(args1, args) # append remaining args to the list
        if(trace >= 2)
            cat0(paste.collapse(names(args)), "\n")
    }
    trace2(trace, "\n")
    args
}
# paste.drop("prefix", "",         drop) returns "prefix,DROP"
# paste.drop("prefix", "x",        drop) returns "prefix,x,DROP,"
# paste.drop("prefix", "x,y",      drop) returns "prefix,x,y,DROP,"
# paste.drop("prefix", c("x","y"), drop) returns "prefix,x,y,DROP,"

paste.drop <- function(prefix, s, drop) {
    s <- paste(s, collapse=",")
    if(nzchar(s))
        paste0(prefix, ",", s, ",", drop)
    else
        paste0(prefix, ",", drop)
}
stopif.unamed.dot <- function(dots, higher.caller, ...) # called from deprefix()
{
    which <- which(names(dots) == "")
    if(length(which)) {
        call <- sprintf("\n       %s\n",
            paste0(strwrap(higher.caller, width=getOption("width"), exdent=10),
                   collapse="\n"))
        dot <- dots[[ which[1] ]]
        env <- parent.frame(2)
        arg <- try(eval(dot, envir=env, enclos=env), silent=TRUE)
        if(is.try.err(arg))
            # fallback to weaker error message "(argument ..1 is unnamed)"
            stop0("Unnamed arguments are not allowed here",
                  " (argument ", as.char(dot), " is unnamed)", call)
        else
            stop0("Unnamed arguments are not allowed here",
                  "\n       The argument's value is ", as.char(arg), call)
    }
}
stopifnot.expanded.dotnames.unique <- function(expanded.names, unexpanded.names,
                                               fname, formals, prefix)
{
    duplicated <- which(duplicated(expanded.names))
    if(length(duplicated) == 0)
        return() # no duplicates
    if(is.null(formals))
        stop0("argument '", unexpanded.names[duplicated[1]],
              "' for ", fname, "() is duplicated")
    else {
        # a little processing is needed because we want to report the
        # error using the unexpanded.names, not the expanded names

        # get the index of the duplicated argument's twin

        duplicated <- duplicated[1]
        for(twin in 1:duplicated)
            if(expanded.names[twin] == expanded.names[duplicated])
                break
        stopifnot(twin < duplicated)

        # get the formal argument matched by the duplicated arguments

        match <- charmatch(sub(prefix, "", expanded.names[duplicated]), formals)

        if(anyNA(match))
            # Dot args are duplicated, but don't match any formal arg.  Probably
            # because e.g. force.xlab is specified but force.xlab is also passed
            # in dots to call.dots (an error in the way call.dots is invoked).
            stop0("argument '", unexpanded.names[duplicated[1]], "' for ",
                  fname, "() is duplicated")

        else if(unexpanded.names[twin] == unexpanded.names[duplicated])
            # dot args are identical and they both match the formal
            stop0("argument '", unexpanded.names[duplicated[1]],
                  "' for ", fname, "() is duplicated")

        else # dot args are not identical but both match the formal
            stop0("'", unexpanded.names[twin], "' and '",
                  unexpanded.names[duplicated], "' both match the '",
                  formals[match[1]], "' argument of ", fname, "()")
    }
}
trace.after.dropkeep <- function(trace, dots)
{
    if(trace >= 2)
        printf("after DROP and KEEP  %s\n", paste.collapse(names(dots)))
}
trace.prolog <- function(trace, prefix, drop, keep, dots, higher.caller)
{
    if(trace >= 2) {
        printf.wrap("TRACE %s", higher.caller)
        printf("\nPREFIX %s\n", prefix)
        printf("DROP %s\n",
            if(is.null(drop)) "NULL"
            else               gsub("\\|>", "\n     >", drop))
        printf("KEEP %s\n",
            if(is.null(keep)) "NULL"
            else               gsub("\\|>", "\n     >", keep))
        names <- names(dots)
        names[which(names=="")] <- "UNNAMED"
        printf("input dotnames       %s\n", paste.collapse(names))
    }
}
