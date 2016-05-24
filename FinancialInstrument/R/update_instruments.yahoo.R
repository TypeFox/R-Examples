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
# $Id: update_instruments.yahoo.R 1638 2014-10-08 03:43:16Z gsee $
#
###############################################################################


#' updates instrument metadata with data from yahoo
#' 
#' Adds/updates information in instrument with data downloaded from yahoo
#' 
#' Although these functions are intended to update the metadata of
#' previously defined instruments, \code{update_instruments.TTR} will
#' define the stocks if they do not already exist.
#'
#' \code{update_instruments.TTR} is only to be used on U.S. stocks denominated in USD.
#'
#' @aliases update_instruments.yahoo update_instruments.TTR
#' @param Symbols can be a vector of instrument names, or, can be \sQuote{all}
#' or \sQuote{stocks} or, for update_instruments.TTR, can be NULL in which case
#' all stocks found with \code{stockSymbols} will be defined
#' @param exchange character vector of names of exchanges. Used in \sQuote{TTR}
#' method. Can be \dQuote{AMEX}, \dQuote{NASDAQ}, or \dQuote{NYSE}
#' @param verbose be verbose?
#' @param silent silence warnings?
#' @return called for side-effect
#' @author Garrett See
#' @seealso \code{\link{update_instruments.instrument}}, 
#'   \code{\link{update_instruments.morningstar}},
#'   \code{\link{update_instruments.masterDATA}},
#'   \code{\link[TTR]{stockSymbols}}, \code{\link{stock}}
#' @references Yahoo! Finance \url{finance.yahoo.com} YahooQuote
#' \url{http://dirk.eddelbuettel.com/code/yahooquote.html} 
#' gummy-stuff.org \url{www.gummy-stuff.org/Yahoo-data.htm} 
#' @examples
#' \dontrun{	
#' 	stock('GS',currency('USD'))
#'  update_instruments.yahoo('GS')
#' 	getInstrument('GS')
#'  update_instruments.TTR('GS')
#'  getInstrument('GS')
#' }
#' @export
update_instruments.yahoo <- function(Symbols=c('stocks','all'), verbose=FALSE ) {
    if (is.null(Symbols) || is.na(Symbols) || !hasArg(Symbols)) Symbols <- 'stocks'
    sym.options <- c('all','stocks')
    symkey <- sym.options[pmatch(Symbols,sym.options)]
    symkey <- na.omit(symkey)[1]
    if (!is.na(symkey)) {
	    if (symkey == 'all' || symkey == 'stocks' || is.null(Symbols)){
            if (symkey == 'all') warning('yahoo can only update stocks.')            
            Symbols <- unique(c(ls_stocks(), ls_instruments_by('src','yahoo')))
        } 
    }
    #make sure it's a vector of instrument names
    if (!is.character(Symbols)) {
        if (verbose) cat('No stocks found to update.\n') 
        return(NULL) #stop('Symbols must be a vector of instrument names, or one of "all", "stocks"')    
    }
    yahoo.syms <- Symbols
    for (i in 1:length(Symbols)) {
        tmp_instr <- try(getInstrument(Symbols[i]),silent=TRUE)
        yahoo.syms[i] <- if (!inherits(tmp_instr, 'try-error') 
                            && !is.null(tmp_instr$src) 
                            && any(names(tmp_instr$src) == 'name')) 
                         { tmp_instr$src$name } else Symbols[i]
    }
    yahoo.syms <- paste(yahoo.syms, collapse=";")
	if (is.null(yahoo.syms) || length(yahoo.syms) == 0) 
        stop('error with symbol names; no Symbols supplied?')
    yahooStuff <- getQuote(yahoo.syms, src="yahoo",
					  what=yahooQF(c("Name", 
                        "Stock Exchange",
					    "Market Capitalization",
					    "Average Daily Volume", 
                        "Earnings/Share", 
					    "EPS Estimate Current Year", 
					    "EPS Estimate Next Year", 
            			"Book Value", "EBITDA",	
                        "52-week Range")))  
#    sym.length <- length(unlist(strsplit(Symbols,";")))    	
    #see yahooQF for available whats
    for (i in 1:length(Symbols)) {
        noNA <- function(x) {
            if (x == 'N/A' || is.na(x)) {NULL} else {x}
        }
        instr <- getInstrument(Symbols[i])
		#Only update stocks from yahoo		
		if (inherits(instr,'stock') || any(instr$src == 'yahoo')) {
		    instr$name=noNA(as.character(yahooStuff[i,2]))
			instr$exchange=noNA(as.character(yahooStuff[i,3]))
			instr$market.cap=noNA(yahooStuff[i,4])
			instr$avg.volume=noNA(suppressWarnings(as.numeric(yahooStuff[i,5])))
			instr$EPS=noNA(suppressWarnings(as.numeric(yahooStuff[i,6])))
			instr$EPS.current.year.est = noNA(suppressWarnings(as.numeric(yahooStuff[i,7])))
			instr$EPS.next.year.est = noNA(suppressWarnings(as.numeric(yahooStuff[i,8])))
			instr$book.value=noNA(suppressWarnings(as.numeric(yahooStuff[i,9])))
			instr$EBITDA=noNA(yahooStuff[i,10])
			instr$range.52wk=noNA(yahooStuff[i,11])
	#		instr$IB=twsSTK(as.character(symdesc[i,1]),'SMART'),

            tclass <- unique(c(class(instr),'instrument'))
            class(instr) <- tclass        
            db <- instr$defined.by
		    if (!is.null(db)) {
		        db <- unlist(strsplit(db,";"))
		        db <- rev(unique(c("yahoo", rev(db))))
		        db <- paste(db,collapse=";") 
		    } else db <- "yahoo"
			instr$defined.by=db 
		    instr$updated=Sys.time()
            
		    assign(Symbols[i], instr, pos=.instrument)
		}
    }        
    Symbols
}

#' @export
#' @rdname update_instruments.yahoo
update_instruments.TTR <- function(Symbols = c("stocks", "all"), 
                                   exchange=c("AMEX","NASDAQ","NYSE"), 
                                   silent=FALSE) {
    if (!suppressWarnings(is.currency.name("USD"))) currency("USD")
    df <- stockSymbols(exchange=exchange)
    if (any(c("stocks", "all") %in% Symbols)) {
        Symbols <- ls_stocks()
    }
    df <- df[df[["Symbol"]] %in% Symbols, ]
    if (nrow(df) == 0) {
        if (!isTRUE(silent)) {
            warning(paste(paste(Symbols,collapse=","), 
                          "not found among those listed on", 
                          paste(exchange,collapse=", ")))
        }
        return(invisible(NULL))
    }
    if (!isTRUE(silent)) {
        cat('defining stocks...\n')
    }
    symout <- NULL    
    for (i in 1:nrow(df)) {
        primary_id <- as.character(df$Symbol[i])
        instr <- try(getInstrument(primary_id, silent = TRUE), silent = TRUE) 
        args <- list()
        arg <- as.list(df[i, ])
        arg$defined.by <- 'TTR'
        if (is.instrument(instr) && !inherits(instr, 'stock')) {
            #make a unique primary_id
            primary_id <- make.names(c(instr$primary_id, 
                                    ls_instruments()),unique=TRUE)[-match(ls_instruments(),
                                        make.names(c(instr$primary_id, ls_instruments()),unique=TRUE))]            
            warning(paste("instrument",instr$primary_id,
                          "is already defined, but not as stock.",
                          "A new instrument", primary_id ,"will be created"))
        } else if (is.instrument(instr)) {
            db <- instr$defined.by
            if (!is.null(db)) {
                db <- unlist(strsplit(db,";"))
                db <- rev(unique(c("TTR", rev(db))))
                db <- paste(db,collapse=";") 
            } else db <- "TTR"
            arg$defined.by <- db
            for (j in 1:length(arg)) {
                instrument_attr(instr$primary_id, names(arg[j]), arg[[j]])
            }
            symout <- c(symout, instr$primary_id)
            next
        }
        arg$primary_id <- primary_id
        arg$currency <- "USD"
        arg$updated <- Sys.time()
        symout <- c(symout, do.call("stock", arg))
    }
    symout 
}


#' Update instrument metadata for ETFs
#'
#' Uses the masterDATA.com list of ETFs and ETNs to update previously defined
#' instruments.
#' 
#' \code{update_instruments.md} is an alias.
#' 
#' MasterDATA classifies each ETF into one of six Fund.Types.  From their
#' website:
#' 
#' US Equity ETF: All constituents trade on a US exchange. Both ProShares and 
#' Rydex sponsor ETFs with the objective of achieving the performance (or a 
#' multiple of the performance) of several major US stock indexes. These ETFs 
#' currently are included in this category despite the fact that their 
#' constituent lists are generally not limited to US stocks.
#' 
#' Global Equity ETF: One or more of the constituents do not trade on a US 
#' Exchange.
#' 
#' Fixed Income ETF:  The constituent list contains government and / or 
#' corporate debt instruments. ETFs with this classification will not be 
#' considered for inclusion in MasterDATA's index / ETF compilation list.
#' 
#' Commodity Based ETF:  This classification of ETF has no constituents but is 
#' structured to reflect the valuation of a commodity such as gold, silver, oil 
#' or interest rates. ETFs with this classification will not be considered for 
#' inclusion in MasterDATA's index / ETF compilation list.
#' 
#' Exchange Traded Notes: A type of unsecured, unsubordinated debt security that 
#' was first issued by Barclays Bank PLC. The purpose of ETNs is to create a 
#' type of security that combines both the aspects of bonds and exchange traded 
#' funds (ETF). Similar to ETFs, ETNs are traded on a major exchange.
#' 
#' @param Symbols character vector of Symbols of ETFs
#' @param silent silence warnings?
#' @return called for side-effect. Each ETF that is updated will be given 
#'   instrument attributes of \dQuote{Name} and \dQuote{Fund.Type}
#' @author Garrett See
#' @references \url{http://masterDATA.com} 
#' (\url{http://www.masterdata.com/helpfiles/ETF_List_Downloads/AllTypes.csv})
#' @seealso \code{\link{update_instruments.yahoo}}, 
#'   \code{\link{update_instruments.instrument}}
#' @examples
#' \dontrun{
#' stock(s <- c("SPY", "DIA"), currency("USD"))
#' update_instruments.masterDATA(s)
#' buildHierarchy(s, "Name", "Fund.Type", "defined.by")
#' }
#' @export
update_instruments.masterDATA <- function(Symbols, silent=FALSE) {
    x <- read.csv("http://www.masterdata.com/helpfiles/ETF_List_Downloads/AllTypes.csv", 
                  stringsAsFactors=FALSE)
    if (missing(Symbols)) Symbols <- unique(c(ls_funds(), ls_stocks()))
    s <- Symbols[Symbols %in% x[["Symbol"]]]
    if (length(s) > 0) {
        # only those that inherit stock or fund
        s <- s[sapply(lapply(s, getInstrument, type=c("stock", "fund"), 
                             silent = TRUE), is.instrument)]
    }
    if (length(s) == 0) {
        if (!isTRUE(silent)) {
            warning("instruments must be defined before this can update them.")
        }
        return(invisible())
    }
    x <- x[!duplicated(x[["Symbol"]]), ]  
    x <- x[x[["Symbol"]] %in% s, ]
    for (i in seq_len(NROW(x))) {
        instrument_attr(x[["Symbol"]][i], "Name", x[["Name"]][i])
        instrument_attr(x[["Symbol"]][i], "Fund.Type", x[["Fund.Type"]][i])
        db <- getInstrument(x[["Symbol"]][i])[["defined.by"]]
        db <- if (is.null(db)) {
            "masterDATA"
        } else paste(unique(c(strsplit(db, ";")[[1]], "masterDATA")), 
                 collapse = ";")
        instrument_attr(x[["Symbol"]][i], "defined.by", db)
        instrument_attr(x[["Symbol"]][i], "updated", Sys.time())
    }
    return(s)
}

#' @export
#' @rdname update_instruments.masterDATA
update_instruments.md <- update_instruments.masterDATA


#' Update instruments with metadata from another instrument.
#'
#' Update instruments with metadata from another instrument.
#' 
#' By default, only attributes that have a value of \code{""} will be given a 
#' new value.
#'
#' If \code{create.new} is \code{TRUE}, then if there are attributes in
#' \code{source_id} that are not in the \code{Symbols}' instrument, those 
#' attributes will be copied to the updated instruments unless they are in 
#' \code{ignore}.
#' 
#' @param Symbols charcter vector of primary_ids or other instrument identifiers.
#' of instruments to be updated.  Alternatively, \code{Symbols} can be an
#' \code{instrument} or list of \code{instrument}s.
#' @param source_id The primary_id (or other identifier) of an instrument, or
#' an instrument.  The \code{source_id} instrument will be used to update the
#' metadata of \code{Symbols}' instruments.
#' @param create.new If FALSE (Default), only attributes that exist but have 
#' empty values will be updated.  If TRUE, new attributes will be created if
#' \code{source_id} has them, but the \code{Symbols} do not.
#' @param ignore vector of names of instrument attributes that should not be
#' copied to the updated instruments.
#' @param assign_i TRUE/FALSE. If TRUE, the updated instruments will be assigned
#' back into the instrument environment.  If FALSE, a list of updated 
#' instruments will be returned
#' @return if \code{isTRUE(assign_i)} a vector of primary_ids of the instruments
#' that were upated.  Otherwise, a list of updated instrument objects.
#' @author Garrett See
#' @seealso \code{\link{update_instruments.yahoo}}, 
#' \code{\link{all.equal.instrument}}
#' @note one way to overwrite attributes of one instrument with those of another
#' is to first set equal to \code{""} those attributes that you want to 
#' overwrite, then use \code{update_instruments.instrument} to copy the 
#' attributes.
#' @examples
#' \dontrun{
#' #rm_instruments()
#' currency("USD")
#' synthetic("SPX", "USD", identifiers=list(yahoo="GSPC"),
#'           tick_size=0.01,
#'          liquidHours="T08:30:00/T15:00:00", 
#'          extraField='something else', 
#'          assign_i=TRUE)
#' stock("SPY", "USD", liquidHours="", assign_i=TRUE)
#' all.equal(getInstrument("SPX"), getInstrument("SPY"))
#' getInstrument("SPY")
#' ## update SPY metadata based on the metadata of SPX
#' ## Only attributes that == "" are updated by default
#' update_instruments.instrument("SPY", "SPX", assign_i=FALSE) #liquidHours
#' update_instruments.instrument("SPY", "SPX", create.new=TRUE,
#'                               ignore=c("identifiers", "type"), 
#'                               assign_i=FALSE)
#' # Although you probably do NOT want to, this will
#' # copy everything new -- including identifiers and type!
#' update_instruments.instrument("SPY", "SPX", create.new=TRUE, ignore=NULL, 
#'                               assign_i=FALSE) 
#' }
#' @export
update_instruments.instrument <- function(Symbols, source_id, create.new=FALSE,
                                          ignore="identifiers", assign_i=TRUE) {
    r <- if (is.instrument(source_id)) { 
        source_id 
    } else getInstrument(source_id)
    if (!is.instrument(r)) {
        stop('source_id is neither an instrument nor the name of an instrument')
    }
    if (is.instrument(Symbols)) {
        Symbols <- list(Symbols)
    }
    out <- lapply(Symbols, function(s) {
        si <- if (!is.instrument(s)) {
            getInstrument(s, silent=TRUE)
        } else s
        if (!is.instrument(si)) {
            warning(paste('could not find instrument"', s, '"Skipping...')) 
            return(NULL)
        }
        all.empty <- do.call(c, lapply(si, function(x) all(x == "")))
        all.empty <- all.empty[!names(all.empty) %in% ignore]
        
        names.empty <- names(all.empty[all.empty])
        for (n in names.empty) {
            if (!is.null(r[[n]])) {
                si[[n]] <- r[[n]]
            }
        }
        if (isTRUE(create.new)) {
            nr <- names(r)
            nr <- nr[!nr %in% ignore]
            nsi <- names(si)
            new.attr <- nr[!nr %in% nsi]
            for (n in new.attr) {
                si[[n]] <- r[[n]]
            }
        }
        db <- si$defined.by
        if (!is.null(db)) {
            db <- unlist(strsplit(db,";"))
            db <- rev(unique(c(r$primary_id, rev(db))))
            db <- paste(db,collapse=";") 
        } else db <- r$primary_id
        si$defined.by <- db 
        si$updated <- Sys.time()
        si
    })
    if (isTRUE(assign_i)) {
        invisible(lapply(out, function(x) {
            if (!is.null(x)) assign(x$primary_id, x, 
                                    pos=.instrument)
        }))
    } else return(out)
    do.call(c, lapply(out, "[[", "primary_id"))
}

