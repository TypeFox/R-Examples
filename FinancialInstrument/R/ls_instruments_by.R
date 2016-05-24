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
# $Id: ls_instruments_by.R 1498 2013-08-25 00:26:39Z gsee $
#
###############################################################################


#' Subset names of instruments
#' 
#' list names of instruments that have an attribute that matches some value
#' 
#' list instruments that have a given attribute level with a given value.
#' 
#' @param what What attribute? (e.g. \dQuote{currency}, \dQuote{type}, 
#'   \dQuote{strike}, etc.)
#' @param value What value must the attribute have? (e.g. \dQuote{EUR},
#' \dQuote{option}, 100, etc.).  If missing or \code{NULL}, the names of all
#'   instruments that have a \code{what} slot will be returned
#' @param pattern only return instruments with \code{pattern} in the name
#' @param match should pattern match names exactly?
#' @param in.slot If the attribute you are looking for is stored inside another
#' slot, this is the name of that slot. (usually "IB")
#' @return chr vector of instrument names
#' @author Garrett See
#' @seealso buildHierarchy, instrument.table, ls_instruments
#' @examples
#' 
#' \dontrun{
#' stock(c("GOOG","INTC"),currency("USD"))
#' synthetic("SnP","USD",src=list(name='^GSPC',src='yahoo'))
#' ls_instruments_by('type','stock')
#' ls_instruments_by("name",NULL,in.slot='src')
#' ls_instruments_by('src',NULL)
#' }
#' @export
ls_instruments_by <- function (what, value, in.slot=NULL, pattern=NULL, match=TRUE) {
    if (length(pattern) > 1 && !match) {
        warning("Using match because length of pattern > 1.")
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
    if (missing(value)) value <- NULL
    tmp_symbols <- NULL 
    for (symbol in symbols) {
        tmp_instr <- try(get(symbol, pos = .instrument),silent=TRUE)
        #TODO: clean this up
        if (is.instrument(tmp_instr)) {
            if (
                if (is.null(value)) { #ls_instruments_by('type',NULL) or ls_instruments_by('name',NULL,'src')
                    if (is.null(in.slot)) { #ls_instruments_by('type',NULL) -- all instruments that have a 'type' element
                        if (!inherits(try(tmp_instr[[what]],silent=TRUE), 
                                      'try-error') && 
                            !is.na(tmp_instr[[what]]) && 
                            !is.null(tmp_instr[[what]])) {TRUE} else {FALSE}
                    } else if (!inherits(try(tmp_instr[[in.slot]][[what]],
                                             silent=TRUE), 'try-error') && 
                        !is.na(tmp_instr[[in.slot]][[what]]) && 
                        !is.null(tmp_instr[[in.slot]][[what]])) {TRUE} else {FALSE}
                } else if (is.null(in.slot)) {
                    if (!is.null(tmp_instr[[what]]) && 
                        !is.na(tmp_instr[[what]]) && 
                        any(tmp_instr[[what]] == value) ) {TRUE} else {FALSE}
                } else { #!is.null(value) && !is.null(in.slot)
                    if (!is.null(tmp_instr[[in.slot]][[what]]) && 
                        !is.na(tmp_instr[[in.slot]][[what]]) && 
                        any(tmp_instr[[in.slot]][[what]] == value)) {TRUE} else {FALSE}
                }
            ) tmp_symbols <- c(tmp_symbols, symbol)
        }    
    }
    tmp_symbols
}

