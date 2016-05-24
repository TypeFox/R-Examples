
ensureDate.internal <- function(x, origin="1970-01-01", dates=NULL) {
    fixme <- factor(x)
    if(all((suppressWarnings(as.numeric(levels(fixme)[fixme]))
            < 1000), na.rm=TRUE)) return(x)
    levels(fixme) <-
        ifelse(grepl("-", levels(fixme)),
               levels(fixme),
               suppressWarnings(format(as.Date(as.integer(levels(fixme)),
                              origin=origin))))
    fixme <- levels(fixme)[fixme]
    return(as.Date(fixme))
}

##' Ensure that x is a Date
##'
##' @aliases ensureDate,Date-method
##' @aliases ensureDate,factor-method
##' @aliases ensureDate,character-method
##' @aliases ensureDate,numeric-method
##' @aliases ensureDate,data.frame-method
##' @aliases ensureDate,list-method
##' @param x Date, fac, chr, num, data.frame, or list that might be or
##' contain something that should be a Date
##' @param origin for Dates that may have lost their attributes,
##' and coerced to integer, the origin for reasserting their
##' Date-ness. default=1970-01-01
##' @param dates for data.frame method, columns of \code{x} that
##' should themselves be ensureDates
##' @return x encoded as a Date \emph{unless all of x is
##' numbers less than 1000,} for datamart numeric survey coding.
setGeneric("ensureDate",
    function(x, origin="1970-01-01", dates=NULL){
        standardGeneric("ensureDate")
    })
setMethod("ensureDate", "Date", ensureDate.internal)
setMethod("ensureDate", "factor", ensureDate.internal)
setMethod("ensureDate", "character", ensureDate.internal)
setMethod("ensureDate", "numeric", ensureDate.internal)
setMethod("ensureDate", "data.frame",
          function(x, origin, dates=c("wave","date")) {
              d <- which(names(x) %in% dates)
              for( i in d )
                  x[,i] <- ensureDate(x[,i])
              return(x)
          })

##' Add a length of NAs before or after x
##' (used in Stacking)
##'
##' @param x a data vector
##' @param before length of NA to insert
##' @param after length of NA to insert
##' @return same vector padded with NAs.
padNA <- function(x, before=0, after=0) {
    if(is.factor(x)) {
        if(any(is.na(x))){
            y <- addNA(factor(x))
            lev <- levels(x)
            if(length(lev)==0)
                lev <- NA
        } else {
            y <- factor(x)
            lev <- levels(y)
        }
        out <- factor(c(rep(NA, before), y, rep(NA, after)))
        out <- factor(out, levels=1:length(lev),
                      labels=lev, exclude=NA)
    } else {
        out <- c(rep(NA, before), x, rep(NA, after))
    }
    return(out)
}
