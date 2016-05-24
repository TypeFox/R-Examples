#' Convert price series to/from notional value
#'
#' \code{Notionalize} multiplies all prices by the contract multiplier
#' \code{Denotionalize} divides all prices by the contract multiplier
#'
#' The mulitplier is only applied to columns with prices.  A column is 
#' considered to be a price column if its name contains \dQuote{Open}, 
#'   \dQuote{High}, \dQuote{Low}, \dQuote{Close}, \dQuote{Bid}, \dQuote{Ask}, 
#'   \dQuote{Trade}, \dQuote{Mid}, or \dQuote{Price} and does not contain
#'   \dQuote{Size}, \dQuote{Sz}, \dQuote{Volume}, \dQuote{Qty}, 
#'   \dQuote{Quantity}, \dQuote{OpInt}, \dQuote{OpenInterest} 
#'   (not case-sensitive)
#' @param x an xts object, or an object that is coercible to xts
#' @param name primary_id of the instrument that has the multiplier;
#'   usually the same as the name of \code{x}
#' @param env environment. where to find \code{x} if only its name is provided
#' @return an object of the same class as \code{x}
#' @author Garrett See
#' @examples
#' \dontrun{
#' source("http://tinyurl.com/download-tblox")
#' getSymbols("CL", src='tblox')
#' define_futures.tblox()
#' tail(Notionalize(CL, "CL"))
#' tail(Denotionalize(Notionalize(CL), "CL"))
#' }
#' @export
#' @rdname Notionalize
Notionalize <- function(x, name, env=.GlobalEnv) {
    if (missing(name)) {
        name <- if (!is.character(x)) {
            deparse(substitute(x))
        } else x
    }
    stopifnot(is.instrument.name(name))
    if (is.character(x)) x <- get(x, pos=env)
    xx <- try.xts(x)
    if (isTRUE(attr(x, "notional"))) {
        # at the end of this function, x is given a "notional" attr with a 
        # value of TRUE. Check here.  If there is a "notional" attr that is 
        # TRUE, then this data has already been notionalized
        warning(paste(name, 'has already been notionalized.', 
            'Set attr(, "notional") to FALSE to allow notionalizing again. \n'))
        return(x)
    }
    # Only add dividends to Open, High, Low, Close, Bid, Ask,
    # Trade, Mid, Price
    cn <- grep("Open|High|Low|Close|Bid|Ask|Trade|Mid|Price", colnames(x), 
               ignore.case=TRUE)    
    sz <- grep("Size|Sz|Volume|Qty|Quantity|OpenInterest|OpInt", colnames(x), 
               ignore.case=TRUE)
    cn <- cn[!cn %in% sz]
    # we'll get a warning if cn is all columns, or no columns; suppress it
    out <- suppressWarnings(cbind(sweep(x[, cn], 1, 
        as.numeric(getInstrument(name)$multiplier), "*"), x[, -cn]))
    # set a flag so that if this function is called 
    # again, it will not try to notionalize something
    # that has already been notionalized.
    attr(x, 'notional') <- TRUE
    reclass(out[, colnames(x)], x) # put back in original order and reclass
}


#' @export
#' @rdname Notionalize
Denotionalize <- function(x, name, env=.GlobalEnv) {
    if (missing(name)) {
        name <- if (!is.character(x)) {
            deparse(substitute(x))
        } else x
    }
    stopifnot(is.instrument.name(name))
    if (is.character(x)) x <- get(x, pos=env)
    xx <- try.xts(x)
    if (!isTRUE(attr(x, "notional"))) {
        warning(paste(name, 'is not notional. Nothing to do.\n'))
        return(x)
    }
    # Only add dividends to Open, High, Low, Close, Bid, Ask,
    # Trade, Mid, Price
    cn <- grep("Open|High|Low|Close|Bid|Ask|Trade|Mid|Price", colnames(x), 
               ignore.case=TRUE)    
    sz <- grep("Size|Sz|Volume|Qty|Quantity|OpenInterest|OpInt", colnames(x), 
               ignore.case=TRUE)
    cn <- cn[!cn %in% sz]
    # we'll get a warning if cn is all columns, or no columns; suppress it
    out <- suppressWarnings(cbind(sweep(x[, cn], 1, 
        as.numeric(getInstrument(name)$multiplier), "/"), x[, -cn]))
    attr(x, 'notional') <- FALSE
    reclass(out[, colnames(x)], x) # put back in original order and reclass    
}

