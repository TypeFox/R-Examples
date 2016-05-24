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
# $Id: ls_expiries.R 1498 2013-08-25 00:26:39Z gsee $
#
###############################################################################


#' show unique expiration dates of instruments
#' 
#' show unique expiration dates of instruments
#' 
#' \code{ls_expires} is an alias. (plural of expires?)
#' 
#' type is currently only implemented for \sQuote{derivative}, \sQuote{future},
#' \sQuote{option}, \sQuote{call} and \sQuote{put} internally, a call is made
#' to the appropriate ls_ function.
#' 
#' @aliases ls_expiries ls_expires
#' @param pattern optional regular expression.
#' @param match exact match?
#' @param underlying_id chr name of underlying or vector of underlying_ids. If
#' NULL, all underlyings will be used
#' @param type chr string name of class that instruments to be returned must
#' inherit.
#' @return named chr vector with length of unique expiration dates of
#' derivatives of class \code{type} and having an underlying_id of
#' \code{underlying_id} if given.
#' @note This should be updated to deal with dates instead of character strings
#' @author Garrett
#' @seealso ls_instruments_by for things like e.g.
#' ls_instruments_by('expires','20110916'), ls_instruments, ls_derivatives,
#' ls_options, ls_calls, buildHierarchy, instrument.table
#' @examples
#' 
#' \dontrun{
#' option_series.yahoo('SPY')
#' option_series.yahoo('DIA',NULL)
#' ls_expiries()
#' 
#' }
#' @export
#' @rdname ls_expiries
ls_expiries <- function(pattern=NULL, match=TRUE, underlying_id=NULL, type='derivative') {
    #if (!is.null(pattern)) underlying_id <- ls_underlyings    
    if (is.null(underlying_id))
        underlying_id <- ls_underlyings(pattern,match)
    symbols <- do.call(eval(paste('ls_',type,"s",sep="")),args=list(pattern=pattern) ) #symbols == all derivatives by default
    dates <- NULL   
    underlyings <- NULL
    for (symbol in symbols) { 
        tmp_instr <- try(get(symbol,pos=.instrument),silent=TRUE)
        if (!is.null(tmp_instr$underlying_id) && any(tmp_instr$underlying_id==underlying_id)) { #the underlying_id of this instr mathces one of the one's we're interested in.
        underlying <- tmp_instr$underlying_id            
            if (is.null(tmp_instr$expires)) { #get value for expiry; may be in 'expires' or 'expiry' slot
                if (!is.null(tmp_instr$expiry)) {
                    expiry <- tmp_instr$expiry
                } else expiry <- NULL
            } else expiry <- tmp_instr$expires
        dates <- c(dates, expiry)
        if (!is.null(expiry)) underlyings <- c(underlyings, underlying)                
        }
        #ll <- list(expiry)
        #names(ll) <- underlying
        #dates <- c(dates, ll)
    }
    #cbind(underlyings,dates[-which(duplicated(underlyings))])
    if(!identical(which(duplicated(dates)),integer(0))) {
        expires <- dates[-which(duplicated(dates))]
        names(expires) <- underlyings[-which(duplicated(dates))]    
    } else {
        expires <- dates
        names(expires) <- underlyings
    }
    expires
#    underlying_id <- underlyings[-which(duplicated(dates))]
#    names(underlying_id) <- dates[-which(duplicated(dates))]    
#    data.frame(underlying_id)    
}

#' @export
#' @rdname ls_expiries
ls_expires <- ls_expiries


#ls_instruments_by('expires','20110916')
#ls_expiries(underlying_id=ls_underlyings(ls_calls())) #Nesting
#ls_expiries('SPY')
#ls_expiries(ls_calls())


