##' Compact a ff vector or ffdf data frame
##'
##' Compact takes a ff vector and tries to use the smallest binary data type for this vector.
##' @name compact
##' @aliases .compact compact.ff compact.ffdf
##' @param x \code{ff} or \code{ffdf} object
##' @param use.na \code{logical} if TRUE the resulting ff vector can contain NA, otherwise this is not checked
##' @param ... other parameters, not actually passed but needed to muffle check.
##' @return compact cloned ff vector, or original if no compacting can be done
.compact <- function(x, use.na=TRUE, ...){
   UseMethod(".compact")
}

##' @rdname compact
##' @importFrom ff vmode clone
##' @S3method .compact ff
.compact.ff <- function(x, use.na=TRUE, ...){
    vm <- which(.vmode == vmode(x)) # which of these is x?
    ## just return the true doubles -- this condition is annoying
    ## spss and r both like to make most numeric cols doubles
    if (vm > 9 && !isTRUE(all.equal(as.integer(x[]),x[]))){
        return(x)
    }
    m <- suppressWarnings(find_min_vmode(x))
    if (m < vm){
        clone(x, vmode=.vmode[m], pattern="./")
    } else {
        x
    }
}

##' @rdname compact
##' @S3method .compact ffdf
.compact.ffdf <- function(x, use.na=TRUE, ...){
    ret <- lapply(ff::physical.ffdf(x), .compact.ff, use.na=use.na)
    res <- do.call(ff::ffdf, c(ret, row.names=NULL))
    close(x)
    res
}



##' Upgrade vmode of ff column while stacking
##'
##' Similar to aux function \code{coerce_to_highest_vmode} in ffbase,
##' except this does range/unique values checking rather than simple
##' type checking.
##'
##' @rdname upgradevmode
##' @param x ff
##' @param y ff or vector
##' @param use.na require a vmode that supports NAs? default=TRUE
##' @importFrom ff vmode clone is.ff
.upgradevmode <- function(x, y, use.na=TRUE){
    vm <- which(.vmode == vmode(x))
    ## just return the true doubles -- this condition is annoying
    if (vm > 9 && !isTRUE(all.equal(as.integer(x[]),x[]))){
        return(x)
    }
    idx <- 1:length(.vmode)
    if (ff::is.factor(x)){
        nu <- length(union(levels(x), unique(y[])))
        r <- c(1, max(3,nu))
    } else {
        ry <- range(y, na.rm=TRUE)
        rx <- range(x, na.rm=TRUE)
        r <- range(rx,ry, na.rm=TRUE)
        if(identical(r, c(-Inf, Inf))) {
            rx <- range(ffbase::unique.ff(x), na.rm=TRUE)
            if(is.ff(y)) {
                ry <- range(ffbase::unique.ff(y), na.rm=TRUE)
            } else {
                ry <- range(y)
            }
            r <- range(rx, ry, na.rm=TRUE)
        }
    }
    if(all(is.na(x))) r <- c(1,2)
    m <- (r[1] >= .vmin) & (r[2] <= .vmax)
    if (isTRUE(use.na)){
        m <- m & is.na(.vNA)
    }
    m <- which(m[idx])[1]
    if (is.na(m)) m <- 2 ## use 'logical' if col is 100% NA
    if (m != vm){
        clone(x, vmode=.vmode[m], pattern="./")
    } else {
        x
    }
}
##' Find minimum vmode for non-ff x
##'
##'
##' @param x any R type
##' @param use.na require that the vmode handle NAs
##' @return int index of ff vmode that will hold x.
##' @importFrom ff is.factor .vmin .vmax .vmode .vNA
find_min_vmode <- function(x, use.na=TRUE) {
    ## boolean     1        ushort      8
    ## logical     2        integer     9
    ## quad        3        single     10
    ## nibble      4        double     11
    ## byte        5        complex    12
    ## ubyte       6        raw        13
    ## short       7        character  14

    if (is.factor(x)){
        r <- c(1, max(3,nlevels(x)))
    } else if (!isTRUE(suppressWarnings(all.equal(as.integer(x[]),x[])))){
        return(11)
    } else {
        r <- suppressWarnings(range(x, na.rm=TRUE))
        if(identical(r, c(-Inf, Inf))) {
            r <- range(unique(x), na.rm=TRUE)
        }
    }
    if(all(is.na(x))) r <- c(1,2)
    m <- (r[1] >= .vmin) & (r[2] <= .vmax)
    if (isTRUE(use.na)){
        m <- m & is.na(.vNA)
    }
    m <- which(m[.vmode])[1]
    if (is.na(m)) m <- 2 ## use 'logical' if col is 100% NA
    m
}
