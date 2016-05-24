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
# $Id: ls_instruments.R 1498 2013-08-25 00:26:39Z gsee $
#
###############################################################################


#' List or Remove instrument objects
#' 
#' display the names of or delete instruments, stocks, options, futures,
#' currencies, bonds, funds, spreads, guaranteed_spreads, synthetics,
#' derivatives, or non-derivatives.
#' 
#' ls functions return the names of all the instruments of the class implied by
#' the function name. rm functions remove the instruments of the class implied
#' by the function name
#' 
#' rm_instruments and rm_non_derivatives will not delete currencies unless the
#' keep.currencies argument is FALSE.
#' 
#' For the rm functions, x can be a vector of instrument names, or nothing.  If
#' \code{x} is missing, all instruments of the relevant type will be removed.
#' 
#' It can be useful to nest these functions to get things like futures
#' denominated in USD.
#' 
#' @aliases ls_instruments ls_stocks ls_options ls_option_series ls_futures
#' ls_future_series ls_currencies ls_non_currencies ls_exchange_rates ls_FX
#' ls_bonds ls_funds ls_spreads ls_guaranteed_spreads ls_synthetics
#' ls_derivatives ls_non_derivatives ls_calls ls_puts rm_instruments rm_stocks
#' rm_options rm_option_series rm_futures rm_future_series rm_currencies
#' rm_exchange_rates rm_FX rm_bonds rm_funds rm_spreads rm_synthetics
#' rm_derivatives rm_non_derivatives
#' @param pattern an optional regular expression.  Only names matching
#' \sQuote{pattern} are returned.
#' @param match return only exact matches?
#' @param verbose be verbose?
#' @param include.series should future_series or option_series instruments be
#' included.
#' @param x what to remove. if not supplied all instruments of relevent class
#' will be removed.  For \code{ls_defined.by} x is the string describing how the
#' instrument was defined.
#' @param keep.currencies If TRUE, currencies will not be deleted.
#' @param includeFX should exchange_rates be included in ls_non_currencies
#' results
#' @return ls functions return vector of character strings corresponding to
#' instruments of requested type rm functions are called for side-effect
#' @author Garrett See
#' @seealso ls_instruments_by, ls_by_currency, ls_by_expiry, ls, rm,
#' instrument, stock, future, option, currency, FinancialInstrument::sort_ids
#' @examples
#' 
#' \dontrun{
#' #rm_instruments(keep.currencies=FALSE) #remove everything from .instrument
#' 
#' # First, create some instruments
#' currency(c("USD", "EUR", "JPY"))
#' #stocks
#' stock(c("S", "SE", "SEE", "SPY"), 'USD')
#' synthetic("SPX", "USD", src=list(src='yahoo', name='^GSPC'))
#' #derivatives
#' option('.SPY', 'USD', multiplier=100, underlying_id='SPY')
#' option_series(root_id="SPY", expires='2011-06-18', callput='put', strike=130)
#' option_series(root_id="SPY", expires='2011-09-17', callput='put', strike=130)
#' option_series(root_id="SPY", expires='2011-06-18', callput='call', strike=130)
#' future('ES', 'USD', multiplier=50, expires='2011-09-16', underlying_id="SPX")
#' option('.ES','USD',multiplier=1, expires='2011-06',strike=1350, right='C', underlying_id='ES')
#' 
#' # Now, the examples
#' ls_instruments() #all instruments
#' ls_instruments("SE") #only the one stock
#' ls_instruments("S", match=FALSE) #anything with "S" in name
#' 
#' ls_currencies()
#' ls_stocks() 
#' ls_options() 
#' ls_futures() 
#' ls_derivatives()
#' ls_puts()
#' ls_non_derivatives()
#' #ls_by_expiry('20110618',ls_puts()) #put options that expire on Jun 18th, 2011
#' #ls_puts(ls_by_expiry('20110618')) #same thing
#' 
#' rm_options('SPY_110618C130')
#' rm_futures()
#' ls_instruments()
#' #rm_instruments('EUR') #Incorrect
#' rm_instruments('EUR', keep.currencies=FALSE) #remove the currency
#' rm_currencies('JPY') #or remove currency like this
#' ls_currencies()
#' ls_instruments()
#' 
#' rm_instruments() #remove all but currencies
#' rm_currencies()
#' 
#' option_series.yahoo('DIA')
#' ls_instruments_by('underlying_id','DIA') #underlying_id must exactly match 'DIA'
#' ls_derivatives('DIA',match=FALSE) #primary_ids that contain 'DIA'
#' rm_instruments()
#' }
#' @export
#' @rdname ls_instruments
ls_instruments <- function(pattern=NULL, match=TRUE, verbose=TRUE) {
    if (length(pattern) > 1 && !match) {
        if (verbose)
            warning("Using match=TRUE because length of pattern > 1.")
        #should I use match?
        #or, ignore pattern and return everything?
        #or, do multiple ls calls and return unique
        match <- TRUE    
    }    
    if (!is.null(pattern) && match) {   #there's a pattern and match is TRUE
        symbols <- ls(.instrument, all.names=TRUE)
        symbols <- symbols[match(pattern,symbols)]
    } else if (!match && length(pattern) == 1) { # pattern is length(1) and don't match
        symbols <- ls(.instrument, all.names=TRUE, pattern=pattern)
    } else if (is.null(pattern)) {  #no pattern
        symbols <- ls(.instrument, all.names=TRUE)
    } # else pattern length > 1 & don't match
    
    is.iname <- is.instrument.name(symbols)
    if (!any(is.iname)) return(NULL)
    symbols[is.iname]
}

#' @export
#' @rdname ls_instruments
ls_stocks <- function(pattern=NULL,match=TRUE) {
    symbols <- ls_instruments(pattern,match)
    tmp_symbols <- NULL            
    for (instr in symbols) {
        tmp_instr <- try(get(instr, pos = .instrument),silent=TRUE)
        if (inherits(tmp_instr, 'stock') && inherits(tmp_instr, 'instrument')) {
            tmp_symbols <- c(tmp_symbols,instr)
        }    
    }
    tmp_symbols
}

#' @export
#' @rdname ls_instruments
ls_options <- function(pattern=NULL,match=TRUE, include.series=TRUE) {
    symbols <- ls_instruments(pattern,match)    
    tmp_symbols <- NULL            
    for (instr in symbols) {
        tmp_instr <- try(get(instr, pos = .instrument),silent=TRUE)
        if (inherits(tmp_instr, 'option') && inherits(tmp_instr, 'instrument')) {
            if (!inherits(tmp_instr, 'option_series') || include.series)
                tmp_symbols <- c(tmp_symbols,instr)
        }    
    }
    tmp_symbols
}

#' @export
#' @rdname ls_instruments
ls_option_series <- function(pattern=NULL,match=TRUE) {
    symbols <- ls_instruments(pattern,match)    
    tmp_symbols <- NULL            
    for (instr in symbols) {
        tmp_instr <- try(get(instr, pos = .instrument),silent=TRUE)
        if (inherits(tmp_instr, 'option_series') && inherits(tmp_instr, 'instrument')) {
            tmp_symbols <- c(tmp_symbols,instr)
        }    
    }
    tmp_symbols
}

#' @export
#' @rdname ls_instruments
ls_futures <- function(pattern=NULL,match=TRUE, include.series=TRUE) {
    symbols <- ls_instruments(pattern,match)
    tmp_symbols <- NULL            
    for (instr in symbols) {
        tmp_instr <- try(get(instr, pos = .instrument),silent=TRUE)
        if (inherits(tmp_instr, 'future') && inherits(tmp_instr, 'instrument')) {
            if (!inherits(tmp_instr, 'future_series') || include.series)
                tmp_symbols <- c(tmp_symbols,instr)
        }    
    }
    tmp_symbols
}

#' @export
#' @rdname ls_instruments
ls_future_series <- function(pattern=NULL,match=TRUE) {
    symbols <- ls_instruments(pattern,match)
    tmp_symbols <- NULL            
    for (instr in symbols) {
        tmp_instr <- try(get(instr, pos = .instrument),silent=TRUE)
        if (inherits(tmp_instr, 'future_series') && inherits(tmp_instr, 'instrument')) {
            tmp_symbols <- c(tmp_symbols,instr)
        }    
    }
    tmp_symbols
}

#' @export
#' @rdname ls_instruments
ls_currencies <- function(pattern=NULL, match=TRUE, includeFX=FALSE) {
    symbols <- ls_instruments(pattern=pattern, match=match)
    tmp_symbols <- NULL
    for (instr in symbols) {
        tmp_instr <- try(get(instr, pos = .instrument),
                         silent=TRUE)
        if (inherits(tmp_instr, 'currency') 
            && inherits(tmp_instr, 'instrument')) {
            if (!inherits(tmp_instr, 'exchange_rate') || isTRUE(includeFX)) {
                tmp_symbols <- c(tmp_symbols,instr)
            }
        }
    }
    tmp_symbols
}

#' @export
#' @rdname ls_instruments
ls_non_currencies <- function(pattern=NULL, match=TRUE, includeFX=TRUE) {
    symbols <- ls_instruments(pattern, match)
    tmp_symbols <- NULL            
    for (instr in symbols) {
        tmp_instr <- try(get(instr, pos = .instrument),
                         silent=TRUE)
        if (!inherits(tmp_instr, 'currency') || 
                (inherits(tmp_instr, 'exchange_rate') && includeFX) ) {
            tmp_symbols <- c(tmp_symbols,instr)
        }
    }
    tmp_symbols
}

#' @export
#' @rdname ls_instruments
ls_exchange_rates <- function(pattern=NULL, match=TRUE) {
    symbols <- ls_currencies(pattern=pattern, match=match, includeFX=TRUE)
    tmp_symbols <- NULL            
    for (instr in symbols) {
        tmp_instr <- try(get(instr, pos = .instrument),
                         silent=TRUE)
        if (inherits(tmp_instr, 'exchange_rate') 
            && inherits(tmp_instr, 'instrument')) {
            tmp_symbols <- c(tmp_symbols,instr)
        }    
    }
    tmp_symbols
}

#' @export
#' @rdname ls_instruments
ls_FX <- ls_exchange_rates

#' @export
#' @rdname ls_instruments
ls_bonds <- function(pattern=NULL,match=TRUE) {
    symbols <- ls_instruments(pattern,match)
    tmp_symbols <- NULL            
    for (instr in symbols) {
        tmp_instr <- try(get(instr, pos = .instrument),silent=TRUE)
        if (inherits(tmp_instr, 'bond') && inherits(tmp_instr, 'instrument')) {
            tmp_symbols <- c(tmp_symbols,instr)
        }    
    }
    tmp_symbols
}

#' @export
#' @rdname ls_instruments
ls_funds <- function(pattern=NULL,match=TRUE) {
    symbols <- ls_instruments(pattern,match)
    tmp_symbols <- NULL            
    for (instr in symbols) {
        tmp_instr <- try(get(instr, pos = .instrument),silent=TRUE)
        if (inherits(tmp_instr, 'fund') && inherits(tmp_instr, 'instrument')) {
            tmp_symbols <- c(tmp_symbols,instr)
        }    
    }
    tmp_symbols
}

#' @export
#' @rdname ls_instruments
ls_spreads <- function(pattern=NULL,match=TRUE) {
    symbols <- ls_instruments(pattern,match)
    tmp_symbols <- NULL            
    for (instr in symbols) {
        tmp_instr <- try(get(instr, pos = .instrument),silent=TRUE)
        if (inherits(tmp_instr, 'spread') && inherits(tmp_instr, 'instrument')) {
            tmp_symbols <- c(tmp_symbols,instr)
        }    
    }
    tmp_symbols
}

#' @export
#' @rdname ls_instruments
ls_guaranteed_spreads <- function(pattern=NULL,match=TRUE) {
    symbols <- ls_instruments(pattern,match)
    tmp_symbols <- NULL            
    for (instr in symbols) {
        tmp_instr <- try(get(instr, pos = .instrument),silent=TRUE)
        if (inherits(tmp_instr, 'guaranteed_spread') && inherits(tmp_instr, 'instrument')) {
            tmp_symbols <- c(tmp_symbols,instr)
        }    
    }
    tmp_symbols
}

#' @export
#' @rdname ls_instruments
ls_synthetics <- function(pattern=NULL, match=TRUE) {
    symbols <- ls_instruments(pattern,match)    
    tmp_symbols <- NULL            
    for (instr in symbols) {
        tmp_instr <- try(get(instr, pos = .instrument),silent=TRUE)
        if (inherits(tmp_instr, 'synthetic') && inherits(tmp_instr, 'instrument')) {
            tmp_symbols <- c(tmp_symbols,instr)
        }    
    }
    tmp_symbols
}

#' @export
#' @rdname ls_instruments
ls_ICS <- function(pattern=NULL, match=TRUE) {
    symbols <- ls_instruments(pattern,match)    
    tmp_symbols <- NULL            
    for (instr in symbols) {
        tmp_instr <- try(get(instr, pos = .instrument),silent=TRUE)
        if (inherits(tmp_instr, 'ICS') && inherits(tmp_instr, 'instrument')) {
            tmp_symbols <- c(tmp_symbols,instr)
        }    
    }
    tmp_symbols
}

#' @export
#' @rdname ls_instruments
ls_ICS_roots <- function(pattern=NULL, match=TRUE) {
    symbols <- ls_instruments(pattern,match)
    tmp_symbols <- NULL            
    for (instr in symbols) {
        tmp_instr <- try(get(instr, pos = .instrument),silent=TRUE)
        if (inherits(tmp_instr, 'ICS_root') && inherits(tmp_instr, 'instrument')) {
            tmp_symbols <- c(tmp_symbols,instr)
        }    
    }
    tmp_symbols
}

# should it be ls_yahoo, ls_defined.by.yahoo, or ls_src? something else?
#ls_yahoo <- function(pattern=NULL) {
#instruments defined by yahoo
#    symbols <- ls_instruments(pattern) #TODO: other functions should be updated to get symbols like this too   
#    tmp_symbols <- NULL
#    for (instr in symbols) {
#        tmp_instr <- try(get(instr, pos = .instrument),silent=TRUE)
#        if ( is.instrument(tmp_instr) && !is.null(tmp_instr$defined.by) )  {        
#            dby <- unlist(strsplit( tmp_instr$defined.by,";"))    
#            if (any(dby == "yahoo" )) 
#                tmp_symbols <- c(tmp_symbols, instr)
#        }
#    }
#    tmp_symbols
#}

#ls_IB <- function(pattern=NULL) {
#instruments defined by IB
#    symbols <- ls_instruments(pattern)   
#    tmp_symbols <- NULL
#    for (instr in symbols) {
#        tmp_instr <- try(get(instr, pos = .instrument),silent=TRUE)
#        if ( is.instrument(tmp_instr) && !is.null(tmp_instr$defined.by) )  {        
#            dby <- unlist(strsplit( tmp_instr$defined.by,";"))    
#            if (any(dby == "IB" )) tmp_symbols <- c(tmp_symbols,instr)
#        }
#    }
#    tmp_symbols
#}


#ls_defined.by <- function(x, pattern=NULL) {
#	symbols <- ls_instruments(pattern)
#	tmp_symbols <- NULL
#	for (symbol in symbols) {
#		tmp_instr <- try(get(symbol, pos=.instrument),silent=TRUE)
#		if (is.instrument(tmp_instr) && !is.null(tmp_instr$defined.by) ) {
#			dby <- unlist(strsplit( tmp_instr$defined.by,";"))
#			if (any(dby == x)) tmp_symbols <- c(tmp_symbols,symbol)
#		}
#	}
#	tmp_symbols
#}

#' @export
#' @rdname ls_instruments
ls_derivatives <- function(pattern=NULL,match=TRUE) {
    symbols <- ls_instruments(pattern,match)
    #there is currently no derivative class    
    #but check for it in case someone made one    
    tmp_symbols <- NULL
    for (instr in symbols) {
        tmp_instr <- try(get(instr, pos = .instrument),silent=TRUE)
         if (inherits(tmp_instr, 'derivative') || 
                 inherits(tmp_instr, 'option') ||
                 inherits(tmp_instr, 'future') ) {
             tmp_symbols <- c(tmp_symbols,instr)
         }    
    }
    tmp_symbols
}

#' @export
#' @rdname ls_instruments
ls_non_derivatives <- function(pattern=NULL,match=TRUE) {
    symbols <- ls_instruments(pattern,match)
    #there is currently no derivative class
    #but check for it in case someone made one    
    tmp_symbols <- NULL
    for (instr in symbols) {
        tmp_instr <- try(get(instr, pos = .instrument),silent=TRUE)
         if (!inherits(tmp_instr, 'derivative') && 
                 !inherits(tmp_instr, 'option') &&
                 !inherits(tmp_instr, 'future') ) {
             tmp_symbols <- c(tmp_symbols,instr)
         }    
    }
    tmp_symbols
}


#' @export
#' @rdname ls_instruments
ls_calls <- function(pattern=NULL,match=TRUE) {
    symbols <- ls_options(pattern=pattern,match=match)
	tmp_symbols <- NULL
    for (instr in symbols) {
        tmp_instr <- try(get(instr, pos = .instrument),silent=TRUE)
		if (is.instrument(tmp_instr) && inherits(tmp_instr, 'option')) {
			if (!is.null(tmp_instr$callput)) {
				right <- tmp_instr$callput
			}else if(!is.null(tmp_instr$right)) {
				right <- tmp_instr$right
			} else right <- FALSE
			if (right == "C" || right == "Call" ||
				right == "call" || right == "c") {
					tmp_symbols <- c(tmp_symbols,instr)
			}			
		}
    }
    tmp_symbols
}

#' @export
#' @rdname ls_instruments
ls_puts <- function(pattern=NULL,match=TRUE) {
    symbols <- ls_options(pattern=pattern,match=match)
	tmp_symbols <- NULL
    for (instr in symbols) {
        tmp_instr <- try(get(instr, pos = .instrument),silent=TRUE)
		if (is.instrument(tmp_instr) && inherits(tmp_instr, 'option')) {
			if (!is.null(tmp_instr$callput)) {
				right <- tmp_instr$callput
			}else if(!is.null(tmp_instr$right)) {
				right <- tmp_instr$right
			} else right <- FALSE
			if (right == "P" || right == "Put" ||
				right == "put" || right == "p") {
					tmp_symbols <- c(tmp_symbols,instr)
			}			
		}
    }
    tmp_symbols
}


#TODO: add error checking: check to see if .instrument exists 

#' @export
#' @rdname ls_instruments
rm_instruments <- function(x, keep.currencies=TRUE) {
    if (missing(x)) {
       x <- ls_instruments()       
    } 
    if (keep.currencies && !is.null(x)) {
        if(any(is.na(match(x,ls_currencies())))) { #are any of them not a currency
            if (!all(is.na(match(x,ls_currencies())))) #are some of them a currency
                x <- x[!x %in% ls_currencies()] #then take them out of to-be-removed
        } else stop('Use keep.currencies=FALSE to delete a currency')    
    }

    rm(list=x,pos=.instrument)
}

#' @export
#' @rdname ls_instruments
rm_stocks <- function(x) {
    if (missing(x)) {
        x <- ls_stocks()
    }
    rm(list=x[x %in% ls_stocks()], pos=.instrument)
}

#' @export
#' @rdname ls_instruments
rm_options <- function(x) {
    if (missing(x)) {
        x <- ls_options()
    }
    rm(list=x[x %in% ls_options()], pos=.instrument)
}

#' @export
#' @rdname ls_instruments
rm_option_series <- function(x) {
    if (missing(x)) {
        x <- ls_option_series()
    }
    rm(list=x[x %in% ls_option_series()], pos=.instrument)
}

#' @export
#' @rdname ls_instruments
rm_futures <- function(x) {
    if (missing(x)) {
        x <- ls_futures()
    }
    rm(list=x[x %in% ls_futures()], pos=.instrument)
}

#' @export
#' @rdname ls_instruments
rm_future_series <- function(x) {
    if (missing(x)) {
        x <- ls_future_series()
    }
    rm(list=x[x %in% ls_future_series()], pos=.instrument)
}

#' @export
#' @rdname ls_instruments
rm_currencies <- function(x) {
    if (missing(x)) {
        x <- ls_currencies()
    }
    rm(list=x[x %in% ls_currencies()], pos=.instrument)
}   

#' @export
#' @rdname ls_instruments
rm_exchange_rates <- function(x) {
    if (missing(x)) {
        x <- ls_exchange_rates()
    }
    rm(list=x[x %in% ls_exchange_rates()], pos=.instrument)
}

#' @export
#' @rdname ls_instruments
rm_FX <- rm_exchange_rates

#' @export
#' @rdname ls_instruments
rm_bonds <- function(x) {
    if (missing(x)) {
        x <- ls_bonds()
    }
    rm(list=x[x %in% ls_bonds()], pos=.instrument)
}

#' @export
#' @rdname ls_instruments
rm_funds <- function(x) {
    if (missing(x)) {
        x <- ls_funds()
    }
    rm(list=x[x %in% ls_funds()], pos=.instrument)
}

#' @export
#' @rdname ls_instruments
rm_spreads <- function(x) {
    if (missing(x)) {
        x <- ls_spreads()
    }
    rm(list=x[x %in% ls_spreads()], pos=.instrument)
}

#' @export
#' @rdname ls_instruments
rm_synthetics <- function(x) {
    if (missing(x)) {
        x <- ls_synthetics()
    }
    rm(list=x[x %in% ls_synthetics()],pos=.instrument)
}


#' @export
#' @rdname ls_instruments
rm_derivatives <- function(x) {
    if (missing(x)) {
        x <- ls_derivatives()
    }
    rm(list=x[x %in% ls_derivatives()],pos=.instrument)
}

#' @export
#' @rdname ls_instruments
rm_non_derivatives <- function(x, keep.currencies=TRUE) {
    if (missing(x)) {
        x <- ls_non_derivatives()
    }
    if (keep.currencies && !is.null(x)) {
        if(any(is.na(match(x,ls_currencies())))) { #are any of them not a currency
            if (!all(is.na(match(x,ls_currencies())))) #are some of them a currency
                x <- x[-match(ls_currencies(),x)] #then take them out of to-be-removed
        } else stop('Use keep.currencies=FALSE to delete a currency')    
    }
    rm(list=x[x %in% ls_non_derivatives()],pos=.instrument) 
}


