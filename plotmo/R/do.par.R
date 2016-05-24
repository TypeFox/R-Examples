# do.par.R: functions setting par() and for setting the overall caption

# main1 is not called main else would clash with main passed in dots (which
# we ignore but cause an error message).  Likewise for xlab1 and ylab1.

do.par <- function(..., nfigs, caption, main1, xlab1, ylab1, trace,
                   nlines.in.main=if(is.specified(main1)) nlines(main1) else 1,
                   def.cex.main=1,
                   def.font.main=2, # use 1 for compat with plot.lm
                   def.right.mar=.8)
{
    nrows <- ceiling(sqrt(nfigs))

    # Note that the plain old cex argument is used in plotmo only in par()
    # (but we also query it later using par("cex")).
    # We use plain old cex relative to the cex calculated by nrows (so passing
    # cex=1 to plotmo causes no changes, and cex=.8 always makes things smaller).
    # TODO cex.axis etc. should be treated in the same way
    # TODO consider moving this into the dotargs functions, also extend for cex.axis, cex.main
    plain.old.cex <- dot("cex", DEF=1, ...)
    check.numeric.scalar(plain.old.cex)
    cex <- if(nrows == 1)      1
           else if(nrows == 2) .83
           else if(nrows >= 3) .66
    cex <- plain.old.cex * cex

    # set oma to make space for caption if necessary
    stopifnot.string(caption, allow.empty=TRUE, null.ok=TRUE)
    def.oma <- dot("oma", ...)
    if(!is.specified(def.oma)) {
        def.oma <- par("oma")
        def.oma[3] <- max(def.oma[3], # .333 to limit cex adjustmment
                          2 + (plain.old.cex^.333 * nlines(caption)))
    }
    cex.lab <- dot("cex.lab",
                   # make the labels small if multiple figures
                   DEF=if(def.cex.main < 1) .8 * def.cex.main else 1, ...)

    mgp <- # compact title and axis annotations
        if(cex.lab < .6)      c(1,   0.2,  0)
        else if(cex.lab < .8) c(1,   0.25, 0)
        else                  c(1.5, 0.4,  0)

    # margins are small to pack plots in, but make bigger if xlab
    # or ylab specified (note that xlab or ylab equal to NULL means
    # that we will later auto generate them)
    mar <- c(
        if(is.null(xlab1) || (is.specified(xlab1) && nzchar(xlab1)))
            4 else 3,               # bottom
        if(is.null(ylab1) || (is.specified(ylab1) && nzchar(ylab1)))
            3 else 2,               # left
        1.2 * nlines.in.main,       # top
        def.right.mar)              # right

    if(nrows >= 5) # small margins if lots of figures
        mar <- cex * mar

    trace2(trace, "\n")
    call.dots(graphics::par,
        DROP="*",               # drop everything
        KEEP="PREFIX,PAR.ARGS", # except args matching PREFIX and PAR.ARGS
        TRACE=if(trace >= 2) trace-1 else 0,
        SCALAR=TRUE,
        def.mfrow     = c(nrows, nrows),
        def.mgp       = mgp,            # compact title and axis annotations
        def.tcl       = -.3,            # shorten tick length
        def.font.main = def.font.main,
        def.mar       = mar,
        def.oma       = def.oma,

        def.cex.main  = def.cex.main,   # ignored by most plot funcs so do it here
        def.cex.lab   = cex.lab,
        def.cex.axis  = cex.lab,

        force.cex     = cex, # last, overrides any cex set by any arg above
        ...) # any remaining graphic dot args are also processed
}
# call do.par on any graphics args in dots, and return a list of their
# old values so the caller can use on.exit to restore them
do.par.dots <- function(..., trace=0)
{
    dots <- match.call(expand.dots=FALSE)$...
    if(length(dots) == 0)
        return(NULL)
    oldpar <- args <- list()
    env <- parent.frame()
    for(dotname in PAR.ARGS) if(is.dot(dotname, ...)) {
        arg <- list(par(dotname))
        names(arg) <- dotname
        oldpar <- append(oldpar, arg)
        dot.org <- dot(dotname, ...)
        dot <- try(eval(dot.org, envir=env, enclos=env), silent=TRUE)
        if(is.try.err(dot))
            dot <- dot.org
        # TODO consider moving this into the dotargs functions, also extend for cex.axis, cex.main
        # special handling for cex args: we want cex to be relative
        # to the current setting, so e.g cex=1 causes no change
        if(substr(dotname, 1, 3) == "cex") {
            olddot <- par(dotname)
            dot <- dot[[1]] * olddot
        } else if (!(dotname %in% PAR.VEC) && length(dot) != 1)
            dot <- dot[[1]] # similar to handling of argument "scalar" in eval.dotlist
        arg <- list(dot)
        names(arg) <- dotname
        args <- append(args, arg)
    }
    if(length(args)) {
        if(trace >= 2)
            printf.wrap("\npar(%s)\n", list.as.char(args))
        do.call(par, args)
    }
    oldpar # a list of old values of args that were changed, empty if none
}
check.do.par <- function(do.par, nfigs) # auto do.par if null, check is 0,1, or 2
{
    if(is.null(do.par))
        do.par <- nfigs > 1
    if(is.logical(do.par))
        do.par <- as.numeric(do.par)
    stopifnot(length(do.par) == 1)
    if(!is.numeric(do.par) || (do.par != 0 && do.par != 1 &&do.par != 2))
        stop0("do.par must be 0, 1, or 2")
    do.par
}
auto.caption <- function(caption, resp.name, type,
                         model.call, object.name, my.call)
{
    sresponse <- stype <- smodel <- scaption <- smy.call <- ""
    if(!is.null(caption))
        scaption <- sprintf("%s     ", caption)
    # the test against "y" is because "y" may just be a fabricated
    # name created because the actual name was not available
    if(!is.null(resp.name) && resp.name != "y")
        sresponse <- paste0(resp.name, "     ")
    if(type != "response")
        stype <- paste0("type=", type, "     ")
    if(!is.null(model.call)) {
        smodel <- strip.deparse(model.call)
        smodel <- sub("\\(formula=", "(", smodel) # delete formula=
    } else
        smodel <- paste0("model: ", object.name)
    s <- paste0(scaption, sresponse, stype, smodel)
    smy.call <- process.my.call.for.caption(my.call)
    if(nzchar(smy.call))
        s <- paste0(s, if(nzchar(s)) "\n" else "", smy.call)
    s
}
# Call this only after a plot is on the screen to avoid
# an error message "plot.new has not been called yet"

draw.caption <- function(caption, ...)
{
    if(!is.null(caption) && any(nzchar(caption))) {
        # allow use of dot args for caption specs
        cex  <- dot("caption.cex  cex.caption",  DEF=1, NEW=1, ...)
        font <- dot("caption.font font.caption", DEF=1, NEW=1, ...)
        col  <- dot("caption.col  col.caption",  DEF=1, NEW=1, ...)
        line <- dot("caption.line", DEF=1, ...)
        # trim so caption fits
        # strwidth doesn't have units of device coords so work with usr coords
        # TODO the algorithm below is not quite correct
        caption <- strsplit(caption, "\n")[[1]]
        usr <- par("usr") # xmin, xmax, ymin, ymax
        n <- par("mfrow")[2] # number of figures horizontally across page
        avail <- .7 * n * (usr[2] - usr[1])
        strwidth <- max(strwidth(caption))
        if(strwidth > avail) {
            which <- strwidth(caption) > avail
            max <- max(nchar(caption))
            max.nchar  <- max * avail / strwidth
            if(max.nchar < max) { # TODO should always be FALSE but actually isn't
                caption <- substr(caption, 1, max.nchar)
                caption[which] <- paste0(caption[which], "...")
            }
        }
        caption <- paste(caption, collapse="\n")
        mtext(text=caption, line=line, outer=TRUE,
              cex=cex * par("cex")^.333, col=col, font=font)
    }
    caption
}
get.caption <- function(nfigs, do.par, caption, resp.name, type,
                        model.call, object.name, my.call)
{
    stopifnot.string(caption, null.ok=TRUE, allow.empty=TRUE)
    if(nfigs > 1 && do.par && (is.null(caption) || !is.null(my.call)))
        auto.caption(caption, resp.name, type,
                     model.call, object.name, my.call)
    else
        paste0(if(is.null(caption)) "" else caption,
            if(!is.null(caption) && !is.null(my.call)) "\n" else "",
            if(!is.null(my.call)) "" else process.my.call.for.caption(my.call))
}
process.my.call.for.caption <- function(my.call)
{
    s <- ""
    if(!is.null(my.call)) {
        s <- sub("\\(object=", "(", my.call)            # delete object=
        s <- sub(", trace=[-._$[:alnum:]]+",     "", s) # delete trace=xxx
        s <- sub(", SHOWCALL=[-._$[:alnum:]]+", "", s)  # delete SHOWCALL=xxx
    }
    s   # a string, may be ""
}
