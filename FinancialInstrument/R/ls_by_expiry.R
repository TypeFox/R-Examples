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
# $Id: ls_by_expiry.R 1498 2013-08-25 00:26:39Z gsee $
#
###############################################################################

#' list or remove instruments by expiration date
#' 
#' show names of or remove instruments that expire on a given date
#' 
#' \code{ls_by_expiry} will find instruments that have a field named either
#' \dQuote{expiry} or \dQuote{expires} with a value that matches \code{expiry}. 
#'  
#' @aliases ls_by_expiry rm_by_expiry
#' @param expiry expiration date that should correspond to the \sQuote{expires}
#' field of an instrument
#' @param pattern an optional regular expression.  Only names matching
#' \sQuote{pattern} are returned.
#' @param match exact match of pattern?
#' @param x what to remove
#' @return \code{ls_by_expiry} gives a vector of names of instruments that
#' expire on the given expiry. \code{rm_by_expiry} is called for its
#' side-effect.
#' @author Garrett See
#' @seealso \code{\link{ls_instruments}}, \code{\link{ls_options}}, \code{\link{ls_calls}}, 
#' \code{\link{ls_puts}}, \code{\link{ls_futures}}, \code{\link{ls_derivatives}}
#' @examples
#' 
#' \dontrun{
#' ls_by_expiry('20110917')
#' ls_by_expiry('20110917',ls_options())
#' }
#' @export
#' @rdname ls_by_expiry
ls_by_expiry <- function(expiry, pattern=NULL, match=TRUE) {
    if (length(pattern) > 1 && !match) {
        warning("Using match because length of pattern > 1.")
        #should I use match even though it's TRUE?
        #or, ignore pattern and return everything?
        #or, do multiple ls calls and return unique
        match <- TRUE    
    }    

    if (!is.null(pattern) && match) {   #there's a pattern and match is TRUE
        symbols <- ls_instruments()
        symbols <- symbols[match(pattern,symbols)]
    } else if (!match && length(pattern) == 1) { # pattern is length(1) and match is FALSE
        symbols <- ls_instruments(pattern=pattern)
    } else if (is.null(pattern)) {  #no pattern
        symbols <- ls_instruments()
    } # else pattern length > 1 & don't match

    expiry <- gsub("-", "", expiry)
    tmp_symbols <- NULL            
    for (symbol in symbols) {
        tmp_instr <- try(get(symbol, pos = .instrument),silent=TRUE)
        if (is.instrument(tmp_instr) ) {
            if ((!is.null(tmp_instr$expires) && any(gsub("-", "", tmp_instr$expires) == expiry)) ||
                (!is.null(tmp_instr$expiry) && any(gsub("-", "", tmp_instr$expiry) == expiry)) ) {
				tmp_symbols <- c(tmp_symbols,symbol)
    		}
        }    
    }
    tmp_symbols
}

#' @export
#' @rdname ls_by_expiry
rm_by_expiry <- function(x,expiry) {
    if (missing(x)) {
        x <- ls_by_expiry(expiry)
    } else x <- ls_by_expiry(expiry,pattern=x)
    rm(list=x,pos=.instrument)
}
#rm_by_expiry(ls_options(),'20130119')


#TODO: ls_by_underlying
        #if (it's  %in% ls_derivatives())
