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
# $Id: all.equal.instrument.R 1655 2014-11-23 22:53:26Z gsee $
#
###############################################################################


#' instrument all.equal method
#'
#' This is most useful for seeing the difference between two \code{instrument} 
#' objects.
#'
#' @param char.n If length of a character vector is \code{char.n} or less it 
#' will be treated as a single element. A negative value for \code{char.n} will
#' be treated as if it were positive \code{Inf}.
#' @param collapse Only used if a character vector is of length less than 
#' \code{char.n}.  Unless \code{collapse} is \code{NULL}, it will be used in a 
#' call to \code{\link{paste}}.  If \code{collapse} is \code{NULL}, each element 
#' of the character vector will be compared separately.
#' @author Garrett See
#' @seealso \code{\link{getInstrument}}, \code{\link{instrument.table}},
#' \code{\link{buildHierarchy}}
#' @note ALPHA code. Subject to change
#' @keywords internal utilities
#' @examples
#' \dontrun{
#' currency("USD")
#' stock("SPY", "USD", validExchanges=c("SMART", "ARCA", "BATS", "BEX"))
#' stock("DIA", "USD", validExchanges=c("SMART", "ARCA", "ISLAND"), 
#'      ExtraField="something")
#' 
#' all.equal(getInstrument("SPY"), getInstrument("DIA"))
#' all.equal(getInstrument("SPY"), getInstrument("DIA"), char.n=5)
#' all.equal(getInstrument("SPY"), getInstrument("DIA"), char.n=5, collapse=NULL)
#' 
#' all.equal(getInstrument("DIA"), getInstrument("USD"))
#' }
#' @export
all.equal.instrument <- function (target, current, char.n=2, collapse=";", ...) {
    stopifnot(is.instrument(target))
    if (char.n < 0) char.n <- Inf
    msg <- NULL
    if (mode(target) != mode(current)) {
        msg <- paste("Modes: ", mode(target), ", ", mode(current), sep="")
    }       
    if (length(target) != length(current)) {
        msg <- c(msg, paste("Lengths: ", length(target), ", ", length(current), 
                            sep=""))
    }
    nt <- names(target)
    nc <- names(current)
    if (is.null(nt) && !is.null(nc)) {
        msg <- c(msg, "names for current but not for target")
        #shouldn't happen because instruments are named lists
    } else if (is.null(nc) && !is.null(nt)) {
        msg <- c(msg, "names for target but not for current")    
    } else {
        if (!all(nt %in% nc)) {
            msg <- c(msg, paste("Names in target that are not in current: <",
                                paste(nt[!nt %in% nc], collapse=", "), ">"))
        }
        if (!all(nc %in% nt)) {
            msg <- c(msg, paste("Names in current that are not in target: <",
                                paste(nc[!nc %in% nt], collapse=", "), ">"))
        }
    }
    if (!is.instrument(current)) {
        msg <- c(msg, paste("target is ", class(target)[1L],
                            ", current is ", class(current)[1L], sep="")) 
        return(msg)
        #TODO: maybe more comparisons can be done depending on what 
        #      class(current) is
    }
    # Same class?
    tc <- class(target)
    cc <- class(current)
    # all instruments have the instrument class, so don't need to compare it
    tc <- tc[!tc %in% "instrument"]
    cc <- cc[!cc %in% "instrument"]
    if (!isTRUE(all.equal(tc, cc))) {
        if (is.null(collapse)) {
            out <- NULL
            if (!all(tc %in% cc)) { 
                out <- paste("Classes of target that are not classes of current: <",
                             paste(tc[!tc %in% cc], collapse=", "), ">")
            }
            if (!all(cc %in% tc)) {
                out <- c(out, paste("Classes of current that are not classes of target: <",
                                    paste(cc[!cc %in% tc], collapse=", "), ">"))
            }
            msg <- c(msg, out)
        } else {
            msg <- c(msg, paste("Classes: ", paste(paste(tc, collapse=collapse), 
                          paste(cc, collapse=collapse), sep=", "), sep=""))
        } 
    }
    uniqueNames <- function(target, current) {  
        unique(c(names(target), names(current)))
    }
    do.compare <- function(target, current, i) {
        if (!isTRUE(all.equal(target[[i]], current[[i]]))) {
            ti <- target[[i]]
            ci <- current[[i]]
            if (is.null(ti)) ti <- "NULL"
            if (is.null(ci)) ci <- "NULL"
            if (is.list(ti)) {
                unames <- uniqueNames(ti, ci)
                out <- do.call(c, lapply(unames, function(x) {
                    if (length(ti) == 1 && ti == "NULL") {
                        paste("NULL, <", names(ci), ">")
                    } else if (length(ci) == 1 && ci == "NULL") {
                        paste("<", names(ti), ">, NULL")
                    } else do.compare(ti, ci, x)
                }))
                return(paste(i, out, sep="$"))
            }
            if (is.xts(ti)) {
                ae <- all.equal(ti, ci)
                if (!isTRUE(ae)) return(paste(i, ae, sep=": "))
            }
            if (max(length(ti), length(ci)) > char.n && is.character(ti)) {
                out <- NULL
                if (!all(ti %in% ci)) 
                    out <- paste(i, "in target but not in current: <",
                                paste(ti[!ti %in% ci], collapse=", "), ">")
                if (!all(ci %in% ti))
                    out <- c(out, paste(i, "in current but not in target: <",
                                paste(ci[!ci %in% ti], collapse=", "), ">"))
                return(out)
            }
            if (!is.null(collapse)) {
                out <- if (isTRUE(all.equal(ti, ci, check.attributes=FALSE))) {
                    all.equal(ti, ci)
                } else {
                    paste(paste(ti, collapse=collapse), 
                          paste(ci, collapse=collapse), sep=", ")
                }
                return(paste(i, ": ", out, sep=""))
            }
            out <- paste(ti, ci, sep=", ")
            out <- paste(i, ": ", out, sep="")
            return(out)
        } 
    }
    ntc <- uniqueNames(target, current)
    msg <- c(msg, 
             do.call(c, lapply(ntc, function(x) do.compare(target, current, x))))
    if (is.null(msg)) {
        TRUE
    } else msg
}
