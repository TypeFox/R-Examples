##' Stack
##' Combine data frames with column names that do not completely intersect
##'
##'
##' @param df1 a \code{data.frame} or a list
##' @param df2 a \code{data.frame} or a list
##' @param return.data.frame defaults to TRUE and uses
##' \code{\link[plyr]{quickdf}}; regular data.frame checks are not done.
##' @param dates character names of columns that are or should be class
##' \code{\link[base]{Date}}. These columns will be coerced to Date, handling
##' the odd case where they are quoted chr days-since-origin.
##' @param origin Origin for dates that are chr like \code{"15218"}, which
##' ought to be 2011-09-01. This is not uncommon for dates coming from SPSS.
##' @param mixed.chr.factor when an element is mixed factor/character in
##' \code{df1} and \code{df2}, what should the result type be? Default="factor".
##' @param verbose print extra information along the way about each column.
##' @return a list of data vectors suitable for coercion to a \code{data.frame}.
##' Because \code{\link{as.data.frame}} is extremely costly and
##' memory-intensive, and Stacking often involves many such column-wise
##' combinations, avoid this by returning a list.
##' @examples
##' testdf1 <- data.frame(both.int=1:4,
##'                       expand.factor=c("blue", "yellow"),
##'                       mixed.fac.int=factor(letters[1:4]),
##'                       date=as.Date("1983-11-22"),
##'                       df1only=rep(1:2, each=2),
##'                       mixed.fac.chr=I(c("a","b","NA",NA)))
##' testdf2 <- data.frame(both.int=5:24,
##'                       expand.factor=factor(rep(c(1:4, NA), 4)),
##'                       mixed.fac.int=1:4,
##'                       date=as.Date("1981-09-24"),
##'                       df2only=factor(c("c", "d")),
##'                       mixed.fac.chr=c("a","b",NA,"c"))
##'
##' levels(testdf2$mixed.fac.chr) <- letters #overleveled
##' ## put levels in a different order than sort(levels) would do, but
##' ## don't make it an ordered factor. Result needs to be ordered
##' ## to preserve this.
##' levels(testdf2$expand.factor) <- c("green", "blue", "red", "yellow")
##' Stack(testdf1, testdf2)
##' Stack(testdf2, testdf1)
##' @export
Stack <- function(df1, df2, return.data.frame=TRUE, dates=c("wave","date"),
                  origin="1970-01-01",
                  mixed.chr.factor=c("factor","character"), verbose=FALSE) {

    mixed.chr.factor <- match.arg(mixed.chr.factor, c("factor","character"))
    names.both <- intersect(names(df1), names(df2)) # names in both
    names.1 <- setdiff(names(df1), names(df2))   # names in 1 but not 2
    names.2 <- setdiff(names(df2), names(df1))   # names in 2 but not 1

    cols.both <- new.env()
    lapply(names.both, function(j) {
        x <- get(j, df1)
        y <- get(j, df2)
        factorx <- is.factor(x)
        factory <- is.factor(y)
        if (factorx && factory && identical(levels(x), levels(y))) {
            assign(j, factor(levels(x)[c(x, y)], levels=levels(x),
                labels=levels(x)), envir=cols.both)
        } else if (!factorx & !factory) {
            assign(j, c( x, y ), envir=cols.both)
        } else if (j %in% dates) {
            x <- ensureDate(x, origin=origin)
            y <- ensureDate(y, origin=origin)
            assign(j, c(x,y), envir=cols.both)
        } else {
            if (factorx) x <- orderIfNeeded(x, verbose=verbose)
            if (factory) y <- orderIfNeeded(y, verbose=verbose)

            if (factorx && factory) {
                lx <- levels(x); ly <- levels(y)
                ordered <- is.ordered(x) || is.ordered(y)
                newlevels <- vectorIntegrate(lx, ly)
                if (!vecEqual(newlevels, sort(newlevels), order=TRUE)) {
                    warning(paste(sQuote(j),': coerced to ordered factor',
                              sep=""),
                              call.=FALSE)
                    ordered <- TRUE
                }
                if (length(lx)!=length(newlevels) ||
                   length(ly)!=length(newlevels)) {
                warning(paste(sQuote(j),
                    ': factor levels expanded: check levels/codes!', sep=""),
                    call.=FALSE)
                }
                new1 <- levels(x)[x]
                new2 <- levels(y)[y]
                assign(j, factor(c(new1, new2), levels=newlevels,
                                 labels=newlevels, ordered=ordered),
                       envir=cols.both)
            }
            if (factorx && !factory & is.numeric(y) &
               length(unique(na.omit(y))) <= length(levels(x))) {
                newlevels <- levels(x)
                assign(j, factor(c(unclass(x),y),
                                 labels=newlevels),
                       envir=cols.both)
                warning(paste(sQuote(j), ": mixed factor-integer: assigned ",
                              'levels from ', sQuote("df1"), sep=""),
                        call.=FALSE)
            }
            if (factory &&
                is.numeric(x) &&
                length(unique(na.omit(x))) <= length(levels(y))) {
                newlevels <- levels(y)
                assign(j,
                       factor(c(x, unclass(y)),
                              labels=newlevels),
                       envir=cols.both)
                warning(paste(sQuote(j), ": mixed integer-factor: assigned ",
                              'levels from ', sQuote("df2"),sep=""),
                        call.=FALSE)
            }

            ## Mixed character-factor
            if (mixed.chr.factor=="factor" &&
                is.factor(x) & is.character(y)) {
                x <- factor(x) # reduce num levels to those present
                newlevels <- unique(c(levels(x),
                               setdiff(levels(x), unique(y))) )
                new1 <- factor(x, levels=newlevels, labels=newlevels)
                new2 <- factor(y, levels=newlevels, labels=newlevels)
                ## explicitly deal with NA level in order to c() them
                if(any(is.na(x)) | any(is.na(y))){
                    new1 <- addNA(new1)
                    new2 <- addNA(new2)
                    newlevels <- levels(new1)
                }
                if (length(newlevels) > length(levels((x)))) {
                    warning(paste(sQuote(j),
                                  ": new levels added from character ",
                                  "values in ",sQuote("df2"),sep=""),
                            call.=FALSE)
                }
                assign(j, factor(c(new1, new2),
                                 labels=newlevels), envir=cols.both)
            }
            if (mixed.chr.factor=="factor" &&
                is.character(x) & is.factor(y)) {
                y <- factor(y) # reduce number of levels to those present
                newlevels <- unique(c(levels(y),
                               setdiff(levels(y), unique(x))) )

                new1 <- factor(x, levels=newlevels, labels=newlevels)
                new2 <- factor(y, levels=newlevels, labels=newlevels)
                ## explicitly deal with NA level in order to c() them
                if (any(is.na(x)) | any(is.na(y))) {
                    new1 <- addNA(new1)
                    new2 <- addNA(new2)
                    newlevels <- levels(new1)
                }
                if (length(newlevels) > length(levels((y)))) {
                    warning(paste(sQuote(j),
                                  ": new levels added from character ",
                                  "values in ",sQuote("df1"),sep=""),
                            call.=FALSE)
                }
                assign(j, factor(c(new1, new2),
                                 labels=newlevels), envir=cols.both)
            }
            if (mixed.chr.factor=="character" &&
                is.factor(x) & is.character(y)) {
                assign(j, c(levels(x)[x], y),
                       envir=cols.both)
                warning(paste(sQuote(j),": mixed factor-character: coerced to ",
                              "character", sep=""),
                        call.=FALSE)
            }
            if (mixed.chr.factor=="character" &&
                is.factor(y) & is.character(x)) {
                assign(j, c(x, levels(y)[y]),
                       envir=cols.both)
                warning(paste(sQuote(j),": mixed factor-character: coerced to ",
                              "character", sep=""),
                        call.=FALSE)
            }
        }
    })
    ## done with cols in both
    cols.1 <- lapply(df1[names.1], FUN=padNA,
                     after=length(df2[[1]]))
    cols.2 <- lapply(df2[names.2], FUN=padNA,
                     before=length(df1[[1]]))
    out <- list()

    out <- c(as.list(cols.both)[names.both] , cols.1, cols.2)
    if(return.data.frame) out <- plyr::quickdf(out)
    return(out)
}

orderIfNeeded <- function(f, verbose=FALSE) {
    if (is.factor(f) && !is.ordered(f)) {
        if (!vecEqual(levels(f), sort(levels(f)), order=TRUE)) {
            f <- ordered(f)
            if (verbose) {
                warning(paste("Asserting order:",
                              paste(sQuote(levels(f)),collapse=", ")),
                        call.=FALSE)
            }
        }
    }
    return(f)
}

do.stack <- function (x, return.data.frame=TRUE, dates=c("wave","date"),
                  origin="1970-01-01",
                  mixed.chr.factor=c("factor","character"), verbose=FALSE) {
    if (is.list(x) && !is.data.frame(x) && length(x)) {
        out <- x[[1]]
        if (length(x)>1) for (i in 2:length(x)) {
            out <- Stack(out, x[[i]], return.data.frame=FALSE, dates=dates,
                origin=origin, mixed.chr.factor=mixed.chr.factor,
                verbose=verbose)
        }
        if (return.data.frame) out <- as.data.frame(out)
        return(out)
    } else {
        if (verbose) warning("Nothing to stack")
        return(x)
    }
}
