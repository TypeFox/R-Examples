#' @name asWide
#' @title Convert an object to "wide" or "long" form.
#'
#' @include ten.R
#' @include nc.R
#' @include print.R
#' 
#' @param x An object of class \code{ten} or \code{pred}.
#' @param ... Additional arguments (not implemented).
#' 
#' @return
#' A new \code{data.table} is returned, 
#' with the data in 'wide' or 'long' format.
#'  \cr
#' There is one row for each time point.
#'  \cr
#' For a \code{ten} object generated from a \code{numeric} or \code{Surv} object, 
#'  this has columns:
#'  \item{t}{\bold{t}ime.}
#'  \item{e}{number of \bold{e}vents.}
#'  \item{n}{\bold{n}umber at risk.}
#' If derived from a \code{survfit}, \code{coxph} or \code{formula} object, 
#' there are additional columns for \code{e} and \code{n}
#' for \emph{each} covariate group. 
#' 
#' @note
#' Most methods for \code{ten} objects are designed for the 'long' form.
#'
#' @rdname asWide
#' @export
#'
asWide <- function(x, ...) UseMethod("asWide")
#'
#' @rdname asWide
#' @method asWide ten
#' @aliases asWide.ten
#' @export
#' @examples
#' data("bmt", package="KMsurv")
#' require("survival")
#' t1 <- ten(c1 <- coxph(Surv(t2, d3) ~ z3*z10, data=bmt))
#' asWide(t1)
#' 
asWide.ten <- function(x, ...){
    data.table::setkey(x, t)
    t1 <- data.table::data.table("t" = x[, sort(unique(t))])
    if (attr(x, "abbNames")) {
        na1 <- attr(x, "longNames")[, id]
        abbFn <- identity
    } else {
        na1 <- attr(x, "longNames")[, longName]
        abbFn <- as.integer
    }
    cg1 <- seq.int(attr(x, "ncg"))
    res1 <- lapply(cg1, FUN=function(cg1){
        r1 <- data.table::setkey(x[abbFn(cg)==cg1, ncg, by=t], t)
        r1 <- r1[t1, roll=-Inf]
        data.table::set(r1, i=which(is.na(r1$ncg)), j="ncg", value=0)
        r1[, ncg]
    })
    res1 <- data.table::as.data.table(res1)
    ## names for 'n' and 'e' columns
    nne1 <- outer(c("n_", "e_"), na1, paste, sep="")
    data.table::setnames(res1, nne1[1, ])
    res2 <-  lapply(cg1, FUN=function(cg1){
        r1 <- data.table::setkey(x[abbFn(cg)==cg1, e, by=t], t)
        r1 <- r1[t1]
        data.table::set(r1, i=which(is.na(r1$e)), j="e", value=0)
        r1[, e]
    })
    res1[, (nne1[2, ]) := res2]
    ## make no. at risk (total) per time period
    res1[, "n" := rowSums(.SD), .SDcols = grep("n_", names(res1))]
    ## total events per time period
    res1[, "e" := rowSums(.SD), .SDcols = grep("e_", names(res1))]
    ## now add time
    res1[, "t" := t1]
    data.table::setcolorder(res1,
                            c("t", "n", "e", as.vector(nne1)))
    data.table::setattr(res1, "class", c("ten", class(res1)))
    setAttr(res1,
            shape="wide",
            abbNames=attr(x, "abbNames"),
            longNames=attr(x, "longNames"),
            ncg=attr(x, "ncg"),
            call=attr(x, "call"),
            mm=attr(x, "mm"))
    return(res1)
}
#'
#' @rdname asWide
#' @export
#'
asLong <- function(x, ...) UseMethod("asLong")
#'
#' @rdname asWide
#' @method asLong ten
#' @aliases asLong.ten
#' @export
#' @examples
#' asLong(asWide(t1))
#' stopifnot(asLong(asWide(t1)) == ten(ten(t1)))
#' 
asLong.ten <- function(x, ...){
    stopifnot(inherits(x, "ten"))
    ## add no. censored
    nc(x)
    ## n at risk + no. events
    n_ <- grep("n_", names(x))
    e_ <- grep("e_", names(x))
    ## names of covariate groups
    n1 <- sub("n_", "", grep("n_", names(x), value=TRUE))
    l1 <- vector(mode="list", length=length(n1))
    for (i in seq_along(n1)){
        e_n <- grep(paste0("e_", n1[i]), names(x))
        c_n <- grep(paste0("c_", n1[i]), names(x))
        ## which times have at least one event
        ## or one censored observation
        l1[[i]] <- as.logical(rowSums(x[, .SD, .SDcols=c(e_n, c_n)]))
    }
### new structure here
    res1 <- data.table::data.table(
        t=unlist(lapply(l1, function(i) x[i, t])),
        n=unlist(lapply(l1, function(i) x[i, n])),
        e=unlist(sapply(seq_along(n1),
                        function(i) x[, .SD, .SDcols=e_[i]][l1[[i]], ])),
        cg=as.integer(unlist(mapply(rep, n1, sapply(l1, sum)))),
        ncg=unlist(sapply(seq_along(n1),
                          function(i) x[, .SD, .SDcols=n_[i]][l1[[i]], ])))
    data.table::setkey(res1, t)
    data.table::setattr(res1, "class", c("ten", class(res1)))
    setAttr(res1,
            "shape"="long",
            "abbNames"=attr(x, "abbNames"),
            "longNames"=attr(x, "longNames"),
            "ncg"=attr(x, "ncg"),
            "call"=attr(x, "call"),
            "mm"=attr(x, "mm"))
    return(res1)
}
## for R CMD check
.SD <- cg <- cg_ <- P <- eMP <- NULL
abbNames <- longName <- id <- NULL
n <- e <- ncg <- NULL
