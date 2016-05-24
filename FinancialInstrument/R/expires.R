###############################################################################
# R (http://r-project.org/) Instrument Class Model
#
# Copyright (c) 2009-2012
# Peter Carl, Dirk Eddelbuettel, Jeffrey Ryan, 
# Joshua Ulrich, Brian G. Peterson, and Garrett See
#
# This library is distributed under the terms of the GNU Public License (GPL)
# for full details see the file COPYING
#
# $Id: expires.R 1655 2014-11-23 22:53:26Z gsee $
#
###############################################################################

#' extract the correct expires value from an \code{instrument}
#'
#' Currently, there are methods for \code{instrument}, \code{spread},
#'  \code{character}, and \code{xts}
#'
#' Will return either the last expiration date before a given \code{Date}, or 
#' the first expiration date after a given \code{Date} 
#' (if \code{expired==FALSE}).
#' 
#' If an \code{\link{instrument}} contains a value for expires that does not
#' include a day (e.g. "2012-03"), or if the expires value is estimated from
#' a \code{future_series} primary_id, it will be assumed that the 
#' \code{instrument} expires on the first of the month (i.e. if the expires
#' value of an instrument were "2012-03", or if there were no expires value
#' but the suffix_id were "H12", the value returned would be "2012-03-01").
#' Note that most non-energy future_series expire after the first of the month 
#' indicated by their suffix_id and most energy products expire in the month
#' prior to their suffix_id month.
#' 
#' @param x instrument or name of instrument
#' @param ... arguments to be passed to methods
#' @return an expiration \code{Date}
#' @author Garrett See
#' @seealso \code{\link{expires.instrument}}, \code{\link{expires.character}}, 
#'   \code{\link{sort_ids}}
#'   
#' \code{\link{getInstrument}} and \code{\link{buildHierarchy}} to see actual 
#'   values stored in \code{instrument}
#' @examples
#' \dontrun{
#' instr <- instrument("FOO_U1", currency=currency("USD"), multiplier=1,
#'                     expires=c("2001-09-01", "2011-09-01", "2021-09-01"), 
#'                     assign_i=FALSE)
#' #Last value of expires that's not after Sys.Date
#' expires(instr) 
#' # First value of expires that hasn't already passed.
#' expires(instr, expired=FALSE)
#' # last value that's not after 2011-01-01
#' expires(instr, Date="2011-01-01") 
#' # first value that's not before 2011-01-01
#' expires(instr, Date="2011-01-01", expired=FALSE) 
#'
#' ## expires.character
#' expires("FOO_U1") # warning that FOO_U1 is not defined
#' instrument("FOO_U1", currency=currency("USD"), multiplier=1,
#'            expires=c("2001-09-01", "2011-09-01", "2021-09-01"), 
#'            assign_i=TRUE)
#' expires("FOO_U1")
#' }
#' @export
expires <- function(x, ...) {
    UseMethod("expires")
}


#' instrument expires extraction method
#' 
#' Returns either the last expiration date before \code{Date}, or the 
#' first expiration date after \code{Date} (if \code{expired==FALSE}).
#' @param Date Can be a Date or character string.  When \code{expires} is a 
#'   vector, the retuned value will be one of the two values of \code{expires} 
#'   that are closest to \code{Date}. (which one will be determined by the value 
#'   of \code{expired})
#' @param expired TRUE/FALSE. This determines which date will be used when
#'   \code{expires} is a vector.  If \code{expired} is \code{TRUE} the date 
#'   returned will be the last one before \code{Date}.  If \code{expired} is 
#'   \code{FALSE} the first one after \code{Date} will be returned. Note that
#'   if \code{expires} is a single value, \code{expired} will be ignored.
#' @param silent silence warnings?
#' @seealso \code{\link{expires}}
#' @author Garrett See
#' @keywords internal
#' @export
expires.instrument <- function(x, Date, expired=TRUE, silent=FALSE, ...) {
    if (is.instrument(x)) {
        if (missing(Date)) Date <- Sys.Date()
        if (!inherits(Date, "Date")) Date <- as.Date(Date)
        xp <- x[["expires"]]
        if (length(xp) == 0) return(NULL)
        chars <- nchar(as.character(xp))
        if (any(chars %in% c(0:5, 9)) || any(chars > 10)) {
            warning(paste("The following values of 'expires' in instrument", 
                x$primary_id, 
                "could not be converted to Date and will be ignored:", 
                xp[chars %in% c(0:5, 9) | chars > 10]))
            xp <- xp[chars %in% c(6:8, 10)]
        }
        if (any(chars < 8) && !isTRUE(silent)) {
            warning(paste("only Year and Month found", "...", 
                          "assuming expiration occurs on the 1st of the month"))
        }
        dxp <- do.call(c, lapply(xp, function(xx) {
            if (inherits(xx, "Date")) {
                xx
            } else if (nchar(xx) == 10) {
                as.Date(xx)
            } else if (nchar(xx) == 8) {
                as.Date(xx, format="%Y%m%d")
            } else if (nchar(xx) == 7) {
                as.Date(paste(xx, "01", sep="-"))
            } else if (nchar(xx) == 6) {
                as.Date(paste(xx, "01", sep=""), format="%Y%m%d")
            }
        }))
        if (length(dxp) == 1) return(dxp)
        if (isTRUE(expired)) {
            return(last(dxp[dxp <= Date]))
        } else return(first(dxp[dxp >= Date]))
    } else NextMethod("expires")
}


#' character expires extraction method
#' 
#' if no \code{instrument} can be found by the id of \code{x}, or if the 
#' \code{instrument} does not have an \code{expires} attribute, an attempt
#' will be made to infer the year and month of expiration using \code{parse_id}
#' in which case the returned value will be a string of the format 
#' \dQuote{YYYY-MM}.  Presently, \code{Date} and \code{expired} will be ignored 
#' if \code{x} is not the name of an instrument
#' @param Date Can be a Date or character string.  When \code{expires} is a 
#'   vector, the retuned value will be one of the two values of \code{expires} 
#'   that are closest to \code{Date}. (which one will be determined by the value 
#'   of \code{expired}).  
#' @param expired TRUE/FALSE. This determines which date will be used when
#'   \code{expires} is a vector.  If \code{expired} is \code{TRUE} the date 
#'   returned will be the last one before \code{Date}.  If \code{expired} is 
#'   \code{FALSE} the first one after \code{Date} will be returned.
#' @param silent silence warnings?
#' @seealso \code{\link{expires.instrument}}
#' @author Garrett See
#' @keywords internal
#' @export
expires.character <- function(x, Date, expired=TRUE, silent=FALSE, ...) {
    xi <- getInstrument(x, silent=TRUE)
    if (is.instrument(xi)) {
        expires.instrument(xi, Date=Date, expired=expired, silent=silent, 
                           ...=...)
    } else {
        if (!isTRUE(silent)) {
            warning(paste(x, "is not defined ... assuming expiration occurs on",
                          "the 1st of the month implied by suffix_id"))
        }
        pid <- parse_id(x)
        mth <- grep(pid$month, month.abb, ignore.case=TRUE)
        mth <- sprintf("%02d", mth, sep="-")
        as.Date(paste(pid$year, mth, "01", sep="-"))
    }
}


#' spread expires extraction method
#' 
#' \code{x$expires} will be returned if it is not \code{NULL}.  Otherwise, the 
#' (character representation of the) exiration date of the first-to-expire of 
#' the \code{members} will be returned.
#' 
#' @param Date Can be a Date or character string.  When \code{expires} is a 
#'   vector, the retuned value will be one of the two values of \code{expires} 
#'   that are closest to \code{Date}. (which one will be determined by the value 
#'   of \code{expired}).  
#' @param expired TRUE/FALSE. This determines which date will be used when
#'   \code{expires} is a vector.  If \code{expired} is \code{TRUE} the date 
#'   returned will be the last one before \code{Date}.  If \code{expired} is 
#'   \code{FALSE} the first one after \code{Date} will be returned.
#' @seealso \code{\link{expires.instrument}}
#' @author Garrett See
#' @keywords internal
#' export
expires.spread <- function(x, Date, expired=TRUE, silent=FALSE, ...) {
    if (inherits(x, "spread")) {
        if (!is.null(x$expires)) {
            return(expires.instrument(x, Date=Date, expired=expired, 
                                      silent=silent, ...=...))
        }
        members <- if (!is.null(x$memberlist$members)) {
            x$memberlist$members
        } else if (!is.null(x$members)) {
            x$members
        } else {
            if (!isTRUE(silent)) {
                warning(paste("Cannot determine members of x$primary_id"))
            }
            return(NextMethod("expires"))
        }
        return(expires.character(sort_ids(members)[1]), Date=Date, 
               expired=expired, silent=silent, ...=...)
    } else NextMethod("expires")
}


#' xts expires extraction method
#' 
#' determines (or estimates) the expiration from an xts object by either
#' finding the last row that is not \code{NA} or by passing the name/symbol
#' of the xts object to \code{\link{expires.character}}
#'
#' @param src either \dQuote{data} or \dQuote{instrument}. 
#' @return If \code{src} is \dQuote{data}, the returned value will be the 
#'   index of the last price that is not \code{NA} (price is determined by 
#'   \code{quantmod:::getPrice}.  \code{getPrice} arguments \code{symbol} and 
#'   \code{prefer} can be passed through dots.)
#'   
#'   If \code{src} is \dQuote{instrument} the symbol of the xts object will
#'   be passed to \code{\link{expires.character}}
#' @seealso \code{\link{expires.instrument}}, \code{\link{expires}}
#' @author Garrett See
#' @keywords internal
#' @export
expires.xts <- function(x, Date, expired=TRUE, silent=FALSE, 
                        src=c("data", "instrument"), ...) {
    src <- c("data", "instrument")[pmatch(src, c("data", "instrument"))[[1L]]]
    if (src == "data") {
        end(na.omit(getPrice(x, ...)))
    } else if (src == "instrument") {
        Symbol <- deparse(substitute(x))
        expires.character(Symbol, Date=Date, expired=expired, silent=silent, 
                          ...=...)
    } else NextMethod("expires")
}
