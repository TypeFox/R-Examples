"times"<-
function(x, ...)
    chron(times. = x, ...)

"Axis.times" <- 
function(x = NULL, at = NULL, ..., side, labels = NULL)
    axis.times(n = side, x = x, labels = labels, ...)

"Math.times" <-
function(x, ...)
{
    if(.Generic == "round")
        return(round_times(x, ...))
    cl <- class(x)
    class(x) <- NULL
    out <- NextMethod(.Generic)
    class(out) <- cl
    out
}

"Ops.times" <-
function(e1, e2)
{
    if(nargs() == 1) {
        ## unary operators
        val <- switch(.Generic,
                      "-" = -1 * e1,
                      "+" = e1,
                      "!" = !as.logical(e1))
        return(val)
    }
    if(is.character(e1))
        e1 <- chron(times. = e1, format = attr(e2, "format"))
    if(is.character(e2))
        e2 <- chron(times. = e2, format = attr(e1, "format"))
    val <- NextMethod(.Generic)
    boolean <- match(.Generic, c("==", "!=", ">", ">=", "<", "<="),
                     nomatch = 0)
    if(boolean) return(as.logical(val))	
    ## make sure the format attribute wasn't dropped by NextMethod
    ## (p.144 blue book) 
    if(is.null(attr(val, "format"))) {
        if(is.null(fmt <- attr(e1, "format")))
            fmt <- attr(e2, "format")
        attr(val, "format") <- fmt
    }
    if(!inherits(val, .Class))
        class(val) <- c(.Class, class(val))
    val
}

"Summary.times" <-
 function(x, ...)
{
    val <- NextMethod(.Generic)
    if(.Generic == "all" || .Generic == "any")
        return(as.logical(val))
    attr(val, "format") <- attr(x, "format")
    class(val) <- class(x)
    val
}

"[.times" <-
function(x, ..., drop = TRUE)
{
    cl <- class(x)
    class(x) <- NULL
    val <- NextMethod("[")
    attr(val, "format") <- attr(x, "format")
    attr(val, "origin") <- attr(x, "origin")
    class(val) <- cl
    val
}

"[<-.times" <-
function(x, ..., value)
{
    if(!as.logical(length(value)))
        return(x)                       # as per p.104 in the blue book
    if(!is.numeric(value) && !is.character(value) && !all(is.na(value)))
        stop("replacement of/with times objects must be with times objects")
    fmt <- attr(x, "format")
    if(!inherits(value, "times"))
        value <- chron(times. = value, format = rev(fmt)[[1]])
    cl <- class(x)                      # ensure that dates objects have
                                        # equal origins
    class(x) <- class(value) <- NULL
    x <- NextMethod(.Generic)
    attr(x, "format") <- fmt
    class(x) <- cl
    x
}

"[[.times" <-
function(x, ..., drop = TRUE)
{
    cl <- class(x)
    class(x) <- NULL
    val <- NextMethod("[[")
    attr(val, "format") <- attr(x, "format")
    attr(val, "origin") <- attr(x, "origin")
    class(val) <- cl
    val
}

"as.character.times" <-
function(x, ...)
    format(x, ...)

as.data.frame.times <- as.data.frame.vector

"axis.times"<-
function(n, x, add = TRUE, labels, simplify = TRUE, ...)
{
    if(!inherits(x, "times"))
        x <- chron(x)
    bad <- is.na(x) | abs(as.vector(x)) == Inf
    rng <- if(n == 1 || n == 3) par("usr")[1:2] else par("usr")[3:4]
    tmp <- c(rng, as.numeric(x[!bad]))
    rng1 <- diff(range(tmp))
    if (rng1 > 1) fctr <- 1
    else if (rng1 > 1/24) fctr <- 24
    else if (rng1 > 1/1440) fctr <- 1440
    else fctr <- 86400
    tmp <- pretty(fctr*tmp)/fctr
    if (simplify) {
        step <- diff(tmp[1:2])
    	simplify <- step >= 1/1440
    	if (inherits(x, "chron") && step >= 1) class(x) <- class(x)[-1]
    }
    
    att <- attributes(x)
    at.x <- structure(tmp[tmp >= rng[1] & tmp <= rng[2]], format = att$
                      format, origin = att$origin, class = att$class)
    if(missing(labels) || (is.logical(labels) && labels)) 
        labels <- format(at.x, simplify = simplify)
    if(add)
        axis(n, at = at.x, labels = labels, ...)
    invisible(list(n = n, at = at.x, labels = labels))
}

"c.times" <-
function(..., recursive = FALSE)
{
    dots <- list(...)
    is.tms <- unlist(lapply(dots, inherits, "times"))
    n <- length(dots)
    fmt <- attr(dots[[(1:n)[is.tms][1]]], "format")
    if(is.null(fmt))
        fmt <- "h:m:s"
    out <- vector("list", length = n)
    for(i in 1:n) {
        x <- dots[[i]]
        if(!all(is.na(x)))
            x <- convert.times(x)
        out[i] <- list(x)
    }
    out <- times(unlist(out, use.names = FALSE), format = fmt)
    out
}

"convert.times"<-
function(times = NULL, format = "h:m:s", length. = 0, ...)
{
    ## convert time in hours, min and secs into fraction of days
    if(is.null(times) || !as.logical(length(times)))
        return(numeric(length = length.))
    if(is.numeric(times))
        return(times)
    if(!is.character(format)) {
        ## format may be a function
        FUN <- switch(mode(format),
                      name = get(as.character(format), mode = "function"),
                      functions = format,
                      stop(paste("unrecognized format mode",
                                 as.character(format))))
        return(FUN(times, ...))
    }
    fmt <- parse.format(format)
    out <- unpaste(times, sep = fmt$sep, fnames = fmt$periods, nfields = 3)
    hh <- mm <- ss <- as.numeric(rep(NA, length(out$h)))
    ok <- !is.na(out$h) & !is.na(out$m) & !is.na(out$s)
    hh[ok] <- as.numeric(out$h[ok])
    mm[ok] <- as.numeric(out$m[ok])
    ss[ok] <- as.numeric(out$s[ok])
    if(all(is.na(hh) | is.na(mm) | is.na(ss)))
        if(any(!is.na(times)))
            stop(paste("format", format, "may be incorrect"))
        else return(rep(NA, length(times)))
    i <- (hh[ok] < 0 | hh[ok] > 23 | mm[ok] < 0 | mm[ok] > 59 |
          ss[ok] < 0 | ss[ok] >= 60)
    bad <- seq_along(hh)[ok][i]
    if(n.bad <- length(bad)) {
        if(n.bad > 10)
            msg <- paste(n.bad, 
                         "time-of-day entries out of range set to NA")
        else msg <- paste("time-of-day entries out of range in positions",
                          paste(bad, collapse = ","), "set to NA")
        warning(msg)
        hh[bad] <- mm[bad] <- ss[bad] <- NA
        ok[bad] <- FALSE
    }
    out <- 3600 * hh + 60 * mm + ss
    out/(24 * 3600)                     # return days and fraction of days
}

"diff.times"<-
function(x, lag = 1, differences = 1, ...)
{
    ## delete references to time-series
    if(lag < 1 | differences < 1)
        stop("Bad value for lag or differences")
    if(lag * differences >= length(x))
        return(x[0])
    r <- x
    s <- 1:lag
    for(i in 1:differences)
        r <- r[ - s] - r[ - (length(r) + 1 - s)]
    r
}

"format.times"<-
function(x, format. = "h:m:s", simplify = FALSE, ...)
{
    if(!as.logical(length(x)))
        return("")
    if(all(is.na(x)))
        return(rep("NA", length = length(x)))
    if(!is.numeric(x))
        stop(paste(deparse(substitute(x)), "must be numeric"))
    att <- attributes(x)
    if(inherits(x, "times")) {
        if(missing(format.))
            format. <- switch(mode(att$format),
                              character = ,
                              list = rev(att$format)[[1]],
                              name = ,
                              "function" = att$format,
                              NULL = format.,
                              stop("invalid output times format"))
        class(x) <- NULL
    }
    if(!is.character(format.)) {
        ## format may be a function or name
        FUN <- switch(mode(format.),
                      "function" = format.,
                      name = eval(format.),
                      stop(paste("unrecognized time format",
                                 deparse(substitute(format.)))))
        return(FUN(unclass(x), ...))
    }
    else format. <- rev(format.)[1]	
    nas <- is.na(x)
    att$class <- att$format <- att$origin <- NULL
    ## <NOTE>
    ## DJ's design is that
    ##   times greater than 1 day  should format like numerics
    ## To change this (e.g., have times(1.5) format as 36:00:00), simply
    ## comment the code below, and make the corresponding change in
    ## print.times().
    days <- abs(floor(x))
    if(any(days[!nas] > 0)) {
        attributes(x) <- att
        return(format(x))
    }
    ## </NOTE>
    sec <- round(24 * 3600 * abs(x))
    hh <- sec %/% 3600
    mm <- (sec - hh * 3600) %/% 60
    ss <- trunc(sec - hh * 3600 - 60 * mm)
    out <- list(h = substring(paste("0", hh, sep = ""), nchar(paste(hh))), 
		m = substring(paste("0", mm, sep = ""), nchar(paste(mm))),
                s = substring(paste("0", ss, sep = ""), nchar(paste(ss))))
    style <- parse.format(format.)
    o <- style$periods
    if(!simplify)
        out <- paste(out[[o[1]]], out[[o[2]]], out[[o[3]]],
                     sep = style$sep)
    else {
        if(simplify == 1) {
            ## no secs
            o <- o[o != "s"]
            out <- paste(out[[o[1]]], out[[o[2]]], sep = style$sep)
        }
        else out <- out$h
    }
    if(any(x[!nas] < 0))
        out <- paste(ifelse(x < 0, "-", " "), out, sep = "")
    out[nas] <- NA
    out[x == Inf] <- "Inf"
    out[x ==  - Inf] <- "-Inf"
    attributes(out) <- att
    out
}

"format<-.times" <-
function(x, ..., value)
{
    ok <- switch(mode(value),
                 character = ,
                 name = ,
                 "function" = ,
                 list = TRUE,
                 FALSE)
    if(!ok)
        stop(paste("invalid format \"", as.character(value), 
                   "\" in format replacement", sep = ""))
    attr(x, "format") <- value
    x
}

"hist.times" <-
function(x, nclass, breaks, plot = TRUE, probability = FALSE, ...,
         xlab = deparse(substitute(x)), simplify = TRUE)
{
    if(!inherits(x, "times"))
        stop(paste(deparse(substitute(x)), "must be of class chron"))
    cl <- class(x)
    x <- as.numeric(x)
    tt <- NextMethod("hist", plot = FALSE)
    dots <- list(...)
    if(plot) {
        old <- par("xaxt", "yaxt")
        on.exit(old)
        x <- tt$breaks
        y <- if(probability) tt$density else tt$counts
        plot.new()
        plot.window(xlim = range(x), ylim = range(y, 0))
        rect(x[-length(x)], 0, x[-1L], y)
        if(any(cl == "dates"))
            lbl <- format(chron(dates. = tt$breaks), simplify = simplify)
        else
            lbl <- format(chron(times. = tt$breaks), simplify = simplify)
        if(is.null(adj <- dots$adj))
            adj <- par("adj")
        if(is.null(cex <- dots$cex))
            cex <- par("cex")
        if(is.null(font <- dots$font))
            font <- par("font")
        if(is.null(las <- dots$las))
            las <- par("las")
        if(is.null(lab <- dots$lab))
            lab <- par("lab")
        if(is.null(mgp <- dots$mgp))
            mgp <- par("mgp")
        if(is.null(tcl <- dots$tcl))
            tcl <- par("tcl")	
        ## do we plot x axis
        if(is.null(axes <- dots$axes))
            axes <- TRUE
        if(is.null(xaxt <- dots$xaxt))
            xaxt <- par("xaxt")
        if(is.null(yaxt <- dots$yaxt))
            yaxt <- par("yaxt")
        if(is.null(horiz <- dots$horiz))
            horiz <- FALSE
        if(axes) {
            if(horiz) {
                if(xaxt != "n")
                    axis(1, adj = adj, labels = TRUE,
                         cex = cex, font = font, las = las, lab = lab,
                         mgp = mgp, tcl = tcl)
            }
            else if(yaxt != "n")
                axis(2, adj = adj, labels = TRUE,
                     cex = cex, font = font, las = las, lab = lab,
                     mgp = mgp, tcl = tcl) 
            axis(horiz + 1, at = tt$breaks, labels = lbl, adj = adj,
                 cex = cex, font = font, las = las, lab = lab, 
                 mgp = mgp, tcl = tcl)
        }
    }
    invisible(tt)
}

"identify.times" <-
function(x, y, ...)
{
    if(inherits(x, "times"))
        x <- as.numeric(x)
    if(!missing(y) && inherits(y, "times"))
        y <- as.numeric(y)
    NextMethod("identify", ...)
}

"is.na.times" <-
function(x, ...)
{
    x <- as.numeric(x)
    NextMethod("is.na")
}

"lines.times" <-
function(x, y, ...)
{
    nas <- is.na(x)
    xtmp <- x <- x[!nas]
    ytmp <- y <- y[!nas]
    o <- order(x)
    x <- as.numeric(x[o])               # as.numeric ensures times are
                                        # computed
    y <- as.numeric(y[o])
    NextMethod("lines", ...)
    invisible(list(x = xtmp, y = ytmp))
}

"mean.times"<-
function(x, trim = 0, weight = rep(1, length(x)), na.ok = TRUE, ...)
{
    if(!missing(weight) && length(weight) != length(x))
        stop(paste("weights must have same length as",
                   deparse(substitute(x))))
    att <- attributes(x)[c("format", "origin", "class")]
    nas <- is.na(x)
    if(!na.ok && any(nas, is.na(weight)))
        return(structure(NA, format = att$format, origin = att$origin, 
                         class = att$class))
    if(na.ok) {
        x <- x[!nas]
        if(!missing(weight))
            weight <- weight[!nas]
    }
    if(trim > 0) {
        if(trim >= 0.5)
            return(median(x))
        n <- length(x)
        i1 <- floor(trim * n) + 1
        i2 <- n - i1 + 1
        i <- sort.list(x, unique(c(i1, i2)))[i1:i2]
        weight <- weight[i]             # lazy eval makes order of
                                        # assignment important!
        x <- x[i]
    }
    if(any(weight < 0))
        stop("weights must be non-negative")
    if(sm <- sum(weight))
        out <- sum(unclass(x) * (weight/sm))
    else out <- rep(0, length(x))
    structure(out, format = att$format, origin = att$origin,
              class = att$class)
}

"plot.times" <-
function(x, y, ...,
         xlab = deparse(substitute(x)), ylab = deparse(substitute(y)),
         simplify)
{
    if(missing(simplify))
        if(is.null(simplify <- getOption("chron.simplify")))
            simplify <- TRUE
    x.times <- inherits(x, "times")	# is x a times?
    if(missing(y)) {
        x <- sort(x)                    # NA's will be ignored
        y <- seq_along(as.vector(x))
        if(missing(ylab))
            ylab <- "Counts"
    }
    y.times <- inherits(y, "times")	# is y a times?
    dots <- list(...)
    if(is.null(axes <- dots$axes)) axes <- TRUE # do we draw axes? 
    ## only xaxt="n" or yaxt="n" requests in ... are honored!
    if(is.null(req.xaxt <- dots$xaxt) || req.xaxt != "n")
        req.xaxt <- "s"
    if(is.null(req.yaxt <- dots$yaxt) || req.yaxt != "n")
        req.yaxt <- "s"
    old <- par("xaxt", "yaxt")
    on.exit(par(old))
    ## trap graphical pars in ... that affect axis() in addition to plot()
    if(is.null(adj <- dots$adj))
        adj <- par("adj")
    if(is.null(cex <- dots$cex.axis))
        cex <- par("cex")
    if(is.null(col <- dots$col.axis))
        col <- par("col")
    if(is.null(font <- dots$font.axis))
        font <- par("font")
    if(is.null(las <- dots$las))
        las <- par("las")
    if(is.null(lab <- dots$lab))
        lab <- par("lab")
    if(is.null(mgp <- dots$mgp))
        mgp <- par("mgp")
    if(is.null(tcl <- dots$tcl)) tcl <- par("tcl")	
    ## for some plot types we need to sort according to x
    if(!is.null(type <- dots$type))
        if(any(type == c("l", "b", "o"))) {
            xlab; ylab                  # force promises
            nas <- is.na(x)
            o <- order(x[!nas])
            x <- x[!nas][o]
            y <- y[!nas][o]
        }
    xx <- unclass(x)
    yy <- unclass(y)
    if(x.times)
        xaxt <- "n"
    else xaxt <- req.xaxt
    if(y.times)
        yaxt <- "n"
    else yaxt <- req.yaxt
    if(!is.null(l <- dots$log)) {
        if(inherits(x, "dates") && any(l == c("x", "xy", "yx")))
            stop("cannot do logarithmic plot of a dates object")
        if(inherits(y, "dates") && any(l == c("y", "xy", "yx")))
            stop("cannot do logarithmic plot of a chron object")
    }
    ## unfortunately we can't use (easily) NextMethod when y is missing!
    plot.default(xx, yy, xlab = xlab, ylab = ylab, ...,
                 xaxt = xaxt, yaxt = yaxt)
    if(axes) {
        if(req.xaxt == "n")
            par(xaxt = "n")
        else if(x.times)
            axis.times(1, x, simplify = simplify, labels = TRUE,
                       adj = adj, col = col, cex = cex, font = font,
                       las = las, lab = lab, mgp = mgp, tcl = tcl)
        if(req.yaxt == "n")
            par(yaxt = "n")
        else if(y.times)
            axis.times(2, y, simplify = simplify, srt = 90, labels
                       = TRUE, adj = adj, col = col, cex = cex,
                       font = font, las = las, lab = lab, mgp = mgp,
                       tcl = tcl)
    }
    invisible(list(x = x, y = y))
}

points.times <- function(x, y, ...)
{
    xtmp <- x
    ytmp <- y
    x <- as.numeric(x)
    y <- as.numeric(y)
    NextMethod("points", ...)
    invisible(list(x = xtmp, y = ytmp))
}

print.times <-
function(x, digits, quote = FALSE, prefix = "", simplify, ...)
{
    if(!as.logical(length(x))) {
        cat("times(0)\n")
        return(invisible(x))
    }
    if(missing(simplify) &&
       is.null(simplify <- getOption("chron.simplify")))
        simplify <- FALSE
    xo <- x
    ## print whole days (no fraction) as regular integers
    if(all(is.na(x)) || any(x[!is.na(x)] >= 1))
        cat("Time in days:\n")
    x <- format.times(x, simplify = simplify)
    NextMethod("print", quote = quote)
    invisible(xo)
}

quantile.times <- function(x, ...)
{
    fmt <- attr(x, "format")
    orig <- attr(x, "origin")
    cl <- class(x)
    x <- unclass(x)
    out <- structure(NextMethod("quantile"), format = fmt, origin = orig, 
                     class = cl)
    out
}

round_times <-
function(x, units = "days", eps = 1e-10, ...)
{
    att <- attributes(x)[c("format", "origin", "class")]
    if(is.character(units)) {
        idx <- pmatch(units, c("days", "hours", "minutes", "seconds"))
        if(!is.na(idx)) {
            values <- c(1, as.numeric(times(c("01:00:00","00:01:00","00:00:01"))))
            units <- values[idx]
        }
    } 
    if(!inherits(units, "times")) { 
        units <- try(times(units))
        if(inherits(units, "try-error")) 
            stop("cannot coerce units to class: times")
    }
    units <- as.numeric(units)
    out <- units * trunc((as.numeric(x) + units / 2 + eps) / units)
    structure(out, format = att$format, origin = att$origin,
              class = att$class)    
}

"summary.times"<-
function(object, digits = 12, ...)
{
    if(!as.logical(length(object)))
        return(object)
    att <- attributes(object)
    class(object) <- NULL
    y <- as.numeric(object)
    z <- unclass(summary.default(y, digits = digits, ...))
    tmp <- structure(z[1:6], format = att$format, origin = att$origin, 
                     class = att$class)
    z[1:6] <- format(tmp)
    class(z) <- "table"
    z
}

## units can be "days", "hours", "minutes", "seconds" or they can
## be of times class or they can be of a class that can be coerced 
## to "times" class
## e.g. trunc(times("12:13:14"), "minutes")         # same
## e.g. trunc(times("12:13:14"), "min")             # same
## e.g. trunc(times("12:13:14"), times("00:01:00")) # same
## e.g. trunc(times("12:13:14"), "00:01:00")        # same
## e.g. trunc(times("12:13:14"), 1/(24*60))         # same
## e.g. trunc(times("12:13:14"), "00:01:30") # truncate to 90 seconds
trunc.times <-
function (x, units = "days", eps = 1e-10, ...)
{
    att <- attributes(x)[c("format", "origin", "class")]
    if(is.character(units)) {
        idx <- pmatch(units, c("days", "hours", "minutes", "seconds"))
        if(!is.na(idx)) {
            values <- c(1, as.numeric(times(c("01:00:00","00:01:00","00:00:01"))))
            units <- values[idx]
        }
    } 
    if(!inherits(units, "times")) { 
        units <- try(times(units))
        if(inherits(units, "try-error")) 
            stop("cannot coerce units to class: times")
    }
    units <- as.numeric(units)
    out <- units * trunc((as.numeric(x) + eps) / units)
    structure(out, format = att$format, origin = att$origin,
              class = att$class)    
}

unique.times <-
function(x, incomparables = FALSE, ...) 
    x[!duplicated(x, incomparables, ...)]

xtfrm.times <-
function(x)
    as.numeric(x)

## chron 'times' objects: only times
## (no dates here because caught by 'dates' method)
pretty.times <-
function(x, ..., simplify = TRUE)
{
   ## call 'chron' method to get absolute times
   ans <- pretty(chron(dates. = x), ...)
   at <- chron(times. = as.vector(ans))
   ## format.times will revert to numeric format if any > 1
   labels <- if(max(abs(at)) <= 1)      # 'x' might exceed 1
       format(at - floor(at), simplify = simplify)
   else
       format(at, simplify = simplify)
   structure(at, labels = labels)
}
