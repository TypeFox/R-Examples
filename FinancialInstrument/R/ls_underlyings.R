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
# $Id: ls_underlyings.R 1498 2013-08-25 00:26:39Z gsee $
#
###############################################################################


#' show names of underlyings
#' 
#' shows names that are stored in the \code{underlying_id} slot of derivative
#' instruments
#' 
#' first calls \code{ls_derivatives}, then looks for unique
#' \code{underlying_id}s. If no derivatives have been defined, nothing will be
#' returned.
#' 
#' @param pattern an optional regular expression.  Only names matching
#' \sQuote{pattern} are returned.
#' @param match require exact match?
#' @return chr vector of names of unique \code{underlying_id}s
#' @author Garrett See
#' @seealso ls_instruments_by, ls_derivatives, ls_options, ls_futures
#' @examples
#' 
#' \dontrun{
#' ls_underlyings()
#' }
#' @export
ls_underlyings <- function(pattern=NULL, match=TRUE) {
    symbols <- ls_derivatives(pattern, match)
    tmp_symbols <- NULL
    for (symbol in symbols) {
        tmp_instr <- try(get(symbol,pos=.instrument),silent=TRUE)
        #if (is.instrument(tmp_instr))  
        if (!is.null(tmp_instr$underlying_id)) 
            tmp_symbols <- c(tmp_symbols,tmp_instr$underlying_id)
    }
    unique(tmp_symbols)    
}

