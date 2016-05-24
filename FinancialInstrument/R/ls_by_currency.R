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
# $Id: ls_by_currency.R 1498 2013-08-25 00:26:39Z gsee $
#
###############################################################################

#' shows or removes instruments of given currency denomination(s)
#' 
#' ls_ functions get names of instruments denominated in a given currency (or
#' currencies) rm_ functions remove instruments of a given currency
#' 
#' 
#' @aliases ls_by_currency rm_by_currency ls_USD ls_AUD ls_GBP ls_CAD ls_EUR
#' ls_JPY ls_CHF ls_HKD ls_SEK ls_NZD
#' @param currency chr vector of names of currency
#' @param pattern an optional regular expression.  Only names matching
#' \sQuote{pattern} are returned.
#' @param match exact match?
#' @param show.currencies include names of currency instruments in the returned
#' names?
#' @param keep.currencies Do not delete currency instruments when deleting
#' multiple instruments.
#' @param x what to remove. chr vector.
#' @return ls_ functions return vector of instrument names rm_ functions return
#' invisible / called for side-effect.
#' @author Garrett See
#' @seealso ls_instruments, ls_currencies, rm_instruments, rm_currencies,
#' twsInstrument, instrument
#' @examples
#' 
#' \dontrun{
#' #First create instruments
#' currency(c('USD','CAD','GBP')
#' stock(c('CM','CNQ'),'CAD')
#' stock(c('BET','BARC'),'GBP')
#' stock(c('SPY','DIA'),'USD')
#' 
#' #now the examples
#' ls_by_currency(c('CAD','GBP'))
#' 
#' ls_USD()
#' ls_CAD()
#' 
#' #2 ways to remove all instruments of a currency
#' rm_instruments(ls_USD()) 
#' #rm_instruments(ls_GBP(),keep.currencies=FALSE)
#' rm_by_currency( ,'CAD') 
#' #rm_by_currency( ,'CAD', keep.currencies=FALSE)
#' }
#' @export
#' @rdname ls_by_currency
ls_by_currency <- function(currency, pattern=NULL, match=TRUE,show.currencies=FALSE) {
    if (length(pattern) > 1 && !match) {
        warning("Using match because length of pattern > 1.")
        #should I use match?
        #or, ignore pattern and return everything?
        #or, do multiple ls calls and return unique
        match <- TRUE    
    }    

	if (!(currency %in% ls_currencies()) ) {
		warning(paste(currency, 'is not a defined currency', sep=" "))
	}

    if (!is.null(pattern) && match) {   #there's a pattern and match is TRUE
        symbols <- ls_instruments()
        symbols <- symbols[match(pattern,symbols)]
    } else if (!match && length(pattern) == 1) { # pattern is length(1) and match is FALSE
        symbols <- ls_instruments(pattern=pattern)
    } else if (is.null(pattern)) {  #no pattern
        symbols <- ls_instruments()
    } # else pattern length > 1 & don't match
        
    tmp_symbols <- NULL            
    for (symbol in symbols) {
        tmp_instr <- try(get(symbol, pos = .instrument),silent=TRUE)
        if (is.instrument(tmp_instr) && 
          tmp_instr$currency == currency ){    
            tmp_symbols <- c(tmp_symbols,symbol)
        }    
    }
    if (show.currencies) {
      tmp_symbols
    } else if (!is.null(tmp_symbols)) {
		ls_non_currencies(tmp_symbols) 
	} else NULL
}

#' @export
#' @rdname ls_by_currency
rm_by_currency <- function(x,currency,keep.currencies=TRUE) {
    sc <- !keep.currencies #make show.currencies==opposite of keep
    if (missing(x)) {
        x <- ls_by_currency(currency,show.currencies=sc)
    } else x <- ls_by_currency(currency,pattern=x,show.currencies=sc)
    rm(list=x,pos=.instrument)
}

#AUD GBP CAD EUR JPY CHF HKD SEK NZD
#' @export
#' @rdname ls_by_currency
ls_USD <- function(pattern=NULL,match=TRUE,show.currencies=FALSE) {
    ls_by_currency('USD',pattern,match,show.currencies) 
}
#' @export
#' @rdname ls_by_currency
ls_AUD <- function(pattern=NULL,match=TRUE,show.currencies=FALSE) {
    ls_by_currency('AUD',pattern,match,show.currencies) 
}
#' @export
#' @rdname ls_by_currency
ls_GBP <- function(pattern=NULL,match=TRUE,show.currencies=FALSE) {
    ls_by_currency('GBP',pattern,match,show.currencies) 
}
#' @export
#' @rdname ls_by_currency
ls_CAD <- function(pattern=NULL,match=TRUE,show.currencies=FALSE) {
    ls_by_currency('CAD',pattern,match,show.currencies) 
}
#' @export
#' @rdname ls_by_currency
ls_EUR <- function(pattern=NULL,match=TRUE,show.currencies=FALSE) {
    ls_by_currency('EUR',pattern,match,show.currencies) 
}
#' @export
#' @rdname ls_by_currency
ls_JPY <- function(pattern=NULL,match=TRUE,show.currencies=FALSE) {
    ls_by_currency('JPY',pattern,match,show.currencies) 
}
#' @export
#' @rdname ls_by_currency
ls_CHF <- function(pattern=NULL,match=TRUE,show.currencies=FALSE) {
    ls_by_currency('CHF',pattern,match,show.currencies) 
}
#' @export
#' @rdname ls_by_currency
ls_HKD <- function(pattern=NULL,match=TRUE,show.currencies=FALSE) {
    ls_by_currency('HKD',pattern,match,show.currencies) 
}
#' @export
#' @rdname ls_by_currency
ls_SEK <- function(pattern=NULL,match=TRUE,show.currencies=FALSE) {
    ls_by_currency('SEK',pattern,match,show.currencies) 
}
#' @export
#' @rdname ls_by_currency
ls_NZD <- function(pattern=NULL,match=TRUE,show.currencies=FALSE) {
    ls_by_currency('NZD',pattern,match,show.currencies) 
}
