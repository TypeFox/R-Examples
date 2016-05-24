##' @importFrom ff ramclass as.ff ff nrow<-
##' @importFrom bit ri
.ffStack <- function(ffdf, df, verbose=FALSE, ...) {
    df1 <- ffdf
    df <- .preparedf(df)
    df2 <- df
    names.both <- intersect(names(df1), names(df2)) # names in both
    names.1 <- setdiff(names(df1), names(df2))   # names in 1 but not 2
    names.2 <- setdiff(names(df2), names(df1))   # names in 2 but not 1
    nx <- nrow(ffdf) ## need for indices to insert y; we extend length next
    ny <- length(df[[1]])
    nrow(ffdf) <- nx + ny
    xi <- ri(1,nx) ## range index of original x
    yi <- ri(nx+1, nx+ny)

    for (j in names.both) {
        if(verbose==TRUE) warning(j)
        x <- ffdf[[j]]
        y <- df[[j]]
        factorx <- ff::is.factor(x)
        factory <- is.factor(y)

        if (!factorx & !factory) {
            ## both are some numeric type
            ffdf[[j]] <- .stackappend(x, y, nx)
        } else if (factorx && factory && identical(levels(x), levels(y))) {
            ffdf[[j]] <- .stackappend(x, y, nx)
        } else {

            if(factorx & !factory & is.numeric(y) &
               length(unique(na.omit(y))) <= length(levels(x))) {
                newlevels <- levels(x)
                y <- factor(levels(x)[y])
                ffdf[[j]] <- .stackappend(x, y, nx)
                warning(paste(sQuote(j), ": mixed factor-integer: assigned ",
                              'levels from ', sQuote("df1"), sep=""),
                        call.=FALSE)
            }
            if (factory==TRUE &&
                ##test from Stack was: is.numeric(x) &&
                is.null(ramclass(x)) &&
                ## unique.ff provided by ffbase
                ## !! need to use xi index here to not select
                ## !! the 0s introduced by extending nrow of ffdf
                length(unique(na.omit(x[xi]))) <= length(levels(y))) {
                newlevels <- levels(y)
                ## 'upgrade' x to a factor
                ## ? can the extend cases be handled here?
                ffdf[[j]] <- as.ff(factor(c(x[xi], y), labels=newlevels))
                warning(paste(sQuote(j), ": mixed integer-factor: assigned ",
                              'levels from ', sQuote("df2"),sep=""),
                        call.=FALSE)
            }
            if (factorx & factory) {
                ffdf[[j]] <- .expandLevels(x[xi], y, j)
            }
        } ## special types
    } ## cols in both

    ## x now has all of its original columns, possibly with unstable
    ## values in y positions. Fill nx+1:ny with NA.
    if (verbose==TRUE) warning(paste(names.1,"in df1 only"))
    ffdf[yi, names.1] <- NA


    for (j in names.2) { ## create the new columns in ffdf, initdata=NA
        if (verbose==TRUE) warning(paste(j, "in df2 only"))
        if(is.factor(df[[j]])) {
            newlevels <- na.omit(levels(df[[j]]))
            ffdf[[j]] <- suppressWarnings(ff(NA, length=nrow(ffdf),
                                             vmode=.vmode[find_min_vmode(df[[j]])],
                                             newlevels)
                                          )
        } else {
            ffdf[[j]] <- ff(NA, length=nrow(ffdf), vmode=.vmode[find_min_vmode(df[[j]])])
        }
        ffdf[yi,j] <- suppressWarnings(df[[j]])
    }
    ffdf
}
.ffffStack <- function(ffdf, df, verbose=FALSE) {
    df1 <- ffdf
    df2 <- df
    names.both <- intersect(names(df1), names(df2)) # names in both
    names.1 <- setdiff(names(df1), names(df2))   # names in 1 but not 2
    names.2 <- setdiff(names(df2), names(df1))   # names in 2 but not 1
    nx <- nrow(ffdf) ## need for indices to insert y; we extend length next
    ny <- nrow(df)
    nrow(ffdf) <- nx + ny
    xi <- ri(1,nx) ## range index of original x
    yi <- ri(nx+1, nx+ny)

    for (j in names.both) {
        x <- ffdf[[j]]
        y <- df[[j]]
        factorx <- ff::is.factor(x)
        factory <- ff::is.factor(y)

        if (!factorx & !factory) {
            ## both are some numeric type
            ffdf[[j]] <- .stackappend(x, y, nx)
        } else if (factorx && factory && identical(levels(x), levels(y))) {
            ffdf[[j]] <- .stackappend(x, y, nx)
        } else {
            if(factorx & !factory &
               length(na.omit(unique(y)[])) <= length(levels(x))) {
                newlevels <- levels(x)
                y <- factor(levels(x)[y[]])
                ffdf[[j]] <- .stackappend(x, y, nx)
                warning(paste(sQuote(j), ": mixed factor-integer: assigned ",
                              'levels from ', sQuote("df1"), sep=""),
                        call.=FALSE)
            }
            if (factory && !factorx  &&
                ## unique.ff provided by ffbase
                ## !! need to use xi index here to not select
                ## !! the 0s introduced by extending nrow of ffdf
                length(na.omit(unique(x[xi])[])) <= length(levels(y))) {
                newlevels <- levels(y)
                ## 'upgrade' x to a factor
                ## ? can the extend cases be handled here?
                ans  <- try(factor(c(x[xi], y[]), labels=newlevels))
                if(!inherits(ans,'try-error')) {
                    ffdf[[j]] <- ff(ans)
                } else {
                    warning(paste(sQuote(j),": levels incompatible"))
                }
                warning(paste(sQuote(j), ": mixed integer-factor: assigned ",
                              'levels from ', sQuote("df2"),sep=""),
                        call.=FALSE)
            }
            if (factorx & factory) {
                ffdf[[j]] <- .expandLevels(x[xi], y, j)
            }
        } ## special types
    } ## cols in both

    ## x now has all of its original columns, possibly with unstable
    ## values in y positions. Fill nx+1:ny with NA.
    ffdf[yi, names.1] <- NA
    for (j in names.2) { ## create the new columns in ffdf, initdata=NA
        if(ff::is.factor(df[[j]])) {
            newlevels <- na.omit(levels(df[[j]]))
            ffdf[[j]] <- suppressWarnings(ff(NA, length=nrow(ffdf),
                                             vmode=vmode(df[[j]]),
                                             levels=newlevels)
                                          )
        } else {
            ffdf[[j]] <- ff(NA, length=nrow(ffdf), vmode="integer")
        }
        ffdf[yi,j] <- df[,j]
    }
    ffdf
}

##' Merge a data.frame into an ffdf data.frame-alike
##'
##' For fast operations on large data frames, we want to turn them into
##' ffdf. This is a special case of \code{\link{Stack}} where the first
##' arugment is already an ffdf, and we are Stacking another data.frame
##' into it, expanding factor levels as needed, and possibly enlarging
##' vmodes of existing ff columns.
##'
##' @rdname ffStack
##' @docType methods
##' @param ffdf an \code{\link[ff]{ffdf}}
##' @param df a \code{\link[base]{data.frame}}
##' @param verbose print extra information about columns as they stack
##' @param ... further arguments
##' @return An ffdf.
##' @import methods
##' @export
setGeneric("ffStack", function (ffdf, df, verbose=FALSE, ...) {
    standardGeneric("ffStack") })
setOldClass("ffdf")
##' @rdname ffStack
##' @aliases ffStack,ffdf,data.frame-method ffStack,ffdf,list-method ffStack,ffdf,ffdf-method
setMethod("ffStack", signature("ffdf", "data.frame"), .ffStack)
setMethod("ffStack", signature("ffdf", "list"), .ffStack)
setMethod("ffStack", signature("ffdf", "ffdf"), .ffffStack)

##' Prepare a data.frame for as.ffdf to not break
##'
##' Presently only converts char to factor.
##' Any other cleaning of things?
##'
##' @rdname preparedf
##' @param df a data.frame or a list
##' @return \code{df} ready for \code{\link[ff]{as.ffdf}}
##' @S3method .preparedf data.frame
##' @S3method .preparedf list
##' @export
.preparedf <- function(df) UseMethod(".preparedf")
.preparedf.data.frame <- function(df){
    chr.cols <- which(sapply(df, is.character))
    if(length(chr.cols)>0) {
        df[,chr.cols] <- lapply(df[,chr.cols,drop=FALSE], function(col) {
            factor(col, exclude=c("NA",NA))
        })
    }
    df
}
.preparedf.list <- function (df) {
    chr.cols <- which(sapply(df, is.character))
    if (length(chr.cols) > 0) {
        for (c in chr.cols) {
            df[[c]] <- factor(df[[c]],
                              exclude = c(NA, "NA"))
        }
    }
    df
}

##' @importFrom ff as.ff clone
##' @importFrom bit as.which chunk
.stackappend <- function(x, y, nx, adjustvmode=TRUE, ...){
    if (is.null(x)){
        if (is.ff(y)){
            return(clone(y, pattern="./"))
        } else {
            return (if (length(y)) ff::as.ff(y))
        }
        nx <- length(x)
    }
    ##TODO check if x and y are compatible
    len <- nx
    to <- length(y)
    if (!to) return(x)

    if (ff::is.factor(x)){
        levels(x) <- ff::appendLevels(levels(x), levels(y))
    }
    ## Upgrade to a higher vmode if needed
    if(adjustvmode==TRUE){
        x <- .upgradevmode(x=x, y=y)
    }
    for (i in chunk(x, from=(nx+1), to=to+nx, ...)){
        i <- i-nx
        if (!is.ff(y)){
            i <- as.which(i)
        }
        x[(i+nx)] <- y[i]
    }
    x
}


.expandLevels <- function(x, y, j) {
    ## get and order the target levels
    newlevels <- vectorIntegrate(na.omit(levels(x)),
                                 levels(y))
    ## supposed to be able to do with with ff level recoding
    ## I could not get it to do this consistently; reverting to
    ## more expensive new factor/new ff. :/

    ##levels(x) <- appendLevels(levels(x), newlevels)
    ##x <- recodeLevels(x, newlevels)
    y <- factor(y[], levels=newlevels)
    ans <- ff::as.ff(factor(c(levels(x)[x[]], levels(y)[y[]]),
                      levels=newlevels, labels=newlevels))
    warning(paste(sQuote(j),
                  ': factor levels expanded: check levels/codes!',
                  sep=""),
            call.=FALSE)
    ans
}
