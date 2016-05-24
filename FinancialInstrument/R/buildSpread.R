#' Construct a price/level series for pre-defined multi-leg spread instrument
#' 
#' Build price series for spreads, butterflies, or other synthetic instruments, 
#' using metadata of a previously defined synthetic instrument.
#'
#' The spread and all legs must be defined instruments.
#'
#' This function can build multileg spreads such as calendars, butterflies, 
#' condors, etc. However, the returned series will be univariate. It does not 
#' return multiple columns (e.g. \sQuote{Bid}, \sQuote{Ask}, \sQuote{Mid}) like 
#' \code{\link{fn_SpreadBuilder}} does.
#'
#' \code{buildBasket} is an alias
#'
#' TODO: allow for multiplier (divisor) that is a vector.
#' @param spread_id The name of the \code{instrument} that contains members and 
#'   memberratio 
#' @param Dates Date range on which to subset.  Also, if a member's data is not 
#'   available via \code{\link{get}} \code{\link[quantmod]{getSymbols}} will be 
#'   called, and the values of the \code{from} and \code{to} arguments will be 
#'   determined using \code{\link[xts]{.parseISO8601}} on \code{Dates}.
#' @param onelot Should the series be divided by the first leg's ratio?
#' @param prefer Price column to use to build structure.
#' @param auto.assign Assign the spread? If FALSE, the xts object will be 
#'   returned.
#' @param env Environment holding data for members as well as where spread data
#'   will be assigned.
#' @return If \code{auto.assign} is FALSE, a univariate xts object. 
#'   Otherwise, the xts object will be assigned to \code{spread_id} and the 
#'   \code{spread_id} will be returned.
#' @seealso 
#' \code{\link{fn_SpreadBuilder}}
#' \code{\link{spread}} for instructions on defining the spread
#' @author Brian Peterson, Garrett See
#' @note this could also be used to build a basket or a strip by using only 
#'   positive values in memberratio
#' @examples
#' \dontrun{
#' currency("USD")
#' stock("SPY","USD",1)
#' stock("DIA","USD",1)
#' getSymbols(c("SPY","DIA")) 
#'
#' spread("SPYDIA", "USD", c("SPY","DIA"),c(1,-1)) #define it.
#' buildSpread('SPYDIA') #build it.
#' head(SPYDIA)
#'
#' }
#' @rdname buildSpread
#' @export
buildSpread <- function(spread_id, Dates = NULL, onelot=TRUE, prefer = NULL, 
                        auto.assign=TRUE, env=.GlobalEnv) {
    spread_instr <- try(getInstrument(spread_id))
    if (inherits(spread_instr, "try-error") | !is.instrument(spread_instr)) {
        stop(paste("Instrument", spread_instr, " not found, please create it first."))
    }
    if (!inherits(spread_instr, "synthetic")) 
        stop(paste("Instrument", spread_id, " is not a synthetic instrument, please use the symbol of a synthetic."))
    #if (!inherits(try(get(spread_id),silent=TRUE), "try-error") && overwrite==FALSE) #Doesn't work..returns vector of FALSE
	#stop(paste(spread_instr,' price series already exists. Try again with overwrite=TRUE if you wish to replace it.')) 

    spread_currency <- spread_instr$currency 
    spread_mult <- as.numeric(spread_instr$multiplier)
    if (is.null(spread_mult) || spread_mult == 0) spread_mult <- 1
    spread_tick <- spread_instr$tick_size
  
    if (!is.null(Dates)) {
      times <- .parseISO8601(Dates)
      from <- times$first.time
      to <- times$last.time
    }
    spreadseries <- NULL
    for (i in seq_len(length(spread_instr$members))) {
        instr <- try(getInstrument(as.character(spread_instr$members[i])))
        if (inherits(instr, "try-error") || !is.instrument(instr)) 
            stop(paste("Instrument", spread_instr$members[i], " not found, please create it first."))
        instr_currency <- instr$currency
        instr_mult <- as.numeric(instr$multiplier)
        instr_ratio <- spread_instr$memberratio[i]
        instr_prices <- try(get(as.character(spread_instr$members[i]),envir=env),silent=TRUE)
        # If we were able to find instr_prices in env, check to make sure there is data between from and to.
        #if we couldn't find it in env or there's no data between from and to, getSymbols
        if (inherits(instr_prices, "try-error") || length(instr_prices) < 2 || (!is.null(Dates) && length(instr_prices[Dates]) == 0)) {
            if (is.null(Dates)) {
                warning(paste(spread_instr$members[i],"not found in 'env', and no Dates supplied. Trying getSymbols defaults.") )
                instr_prices <- getSymbols(as.character(spread_instr$members[i]),auto.assign=FALSE)
                from <- first(index(instr_prices))
                to <- last(index(instr_prices))
            } else {
                warning(paste("Requested data for", spread_instr$members[i], "not found in 'env'. Trying getSymbols."))
                instr_prices <- getSymbols(as.character(spread_instr$members[i]), from = from, to = to, auto.assign=FALSE)
            }
        }
        if (is.null(Dates)) {
            from <- first(index(instr_prices))
            to <- last(index(instr_prices))
        }
        instr_prices <- instr_prices[paste(from,to,sep="::")]
        ##TODO: if length(prefer > 1), use the first value that exists in colnames(instr_prices)
        ##	i.e. treat prefer as an ordered vector of preferences.
        if (is.null(prefer)) { 
          if (is.HLC(instr_prices)) { 
	        pref='Close'
          } else
          if (has.Mid(instr_prices)) {
	        pref='Mid'
          } else
          if (has.Trade(instr_prices)) {
	        pref='Trade'
          } else
          if (has.Price(instr_prices)) {
	        pref='Price'
          } else pref=colnames(instr_prices)[1]
        } else pref=prefer
        if (!is.logical(instr_prices) || ncol(instr_prices) > 1) instr_prices <- getPrice(instr_prices,prefer=pref)
        if (instr$currency != spread_currency) 
            instr_prices <- redenominate(instr_prices,spread_currency,instr$currency)
        if (ncol(instr_prices) > 1) instr_prices <- getPrice(instr_prices, prefer=pref)
        instr_norm <- instr_prices * instr_mult * instr_ratio
        colnames(instr_norm) <- paste(as.character(spread_instr$members[i]), 
            prefer, sep = ".")
        if (is.null(spreadseries)) {
            spreadseries <- instr_norm
            last.from <- start(spreadseries)
            first.to <- end(spreadseries)        
        } else {
            spreadseries = merge(spreadseries, instr_norm)
            last.from <- max(last.from, start(instr_norm)) 
            first.to <- min(first.to, end(instr_norm))
        }    
    }
    spreadseries <- spreadseries[paste(last.from, first.to, sep="/")]
    spreadseries <- na.locf(spreadseries,na.rm=TRUE)
    #browser()
    spreadlevel = xts(rowSums(spreadseries),order.by=index(spreadseries)) #assumes negative memberratio values for shorts in 'memberratio'
    if (onelot) 
        spreadlevel = spreadlevel/abs(spread_instr$memberratio[1]) #abs() takes care of things like a crack spread which is -3:2:1.
    colnames(spreadlevel) <- paste(spread_id,pref,sep='.')
    #Divide by multiplier and round according to tick_size of spread_instr
    if (is.null(spread_tick) || spread_tick == 0) ret <- spreadlevel*spread_mult
    else ret <- round((spreadlevel * spread_mult) / spread_tick, spread_tick) * spread_tick
    if (auto.assign) {
        assign(spread_id, ret, pos=env)
        ret <- spread_id
    }
    ret
}

#' @rdname buildSpread
#' @export
buildBasket <- buildSpread


#' Calculate prices of a spread from 2 instruments.
#'
#' Given 2 products, calculate spread values for as many columns as practicable. 
#'
#' \code{prod1} and \code{prod2} can be the names of instruments, or the xts 
#' objects themselves. Alternatively, \code{prod2} can be omitted, and a vector 
#' of 2 instrument names can be given to \code{prod1}. See the last example for 
#' this usage.
#' 
#' If \code{prod1} and \code{prod2} are names (not xts data), it will try to get 
#' data for \code{prod1} and \code{prod2} from \code{env} (.GlobalEnv by 
#' default).  If it cannot find the data, it will get it with a call to 
#' getSymbols. Prices are multiplied by multipliers and exchange rates to get 
#' notional values in the currency specified.  The second leg's notional values 
#' are multiplied by \code{ratio}.  Then the difference is taken between the 
#' notionals of leg1 and the new values for leg2.
#' 
#' \sQuote{make.index.unique} uses the xts function \code{make.index.unique} 
#' \sQuote{least.liq} subsets the spread time series, by using the timestamps 
#' of the leg that has the fewest rows.
#' \sQuote{duplicated} removes any duplicate indexes.
#' \sQuote{price.change} only return rows where there was a price change in the 
#' Bid, Mid or Ask Price of the spread.
#'
#' @param prod1 chr name of instrument that will be the 1st leg of a 2 leg 
#' spread (Can also be xts data for first product)
#' @param prod2 chr name of instrument that will be the 2nd leg of a 2 leg 
#' spread (Can also be xts data for second product)
#' @param ratio Hedge ratio. Can be a single number, or a vector of same length 
#'   as data.
#' @param currency chr name of currency denomination of the spread
#' @param from from Date to pass through to getSymbols if needed.
#' @param to to Date to pass through to getSymbols if needed.
#' @param session_times ISO-8601 time subset for the session time, in GMT, in 
#'   the format 'T08:00/T14:59'
#' @param notional TRUE/FALSE. Should the prices be multiplied by contract 
#'   multipliers before calculating the spread?
#' @param unique_method method for making the time series unique
#' @param auto.assign If \code{TRUE} (the default) the constructed spread will 
#'   be stored in symbol created with \code{\link{make_spread_id}}. instrument 
#'   metadata will also be created and stored with the same primary_id.
#' @param env If \code{prod1} and \code{prod1} are character, this is where to 
#'   \code{get} the data.  Also, if \code{auto.assign} is \code{TRUE} this is 
#'   the environment in which to store the data (.GlobalEnv by default) 
#' @param silent silence warnings? (FALSE by default)
#' @param \dots other arguments to pass to \code{getSymbols} and/or 
#'   \code{\link{make_spread_id}}
#' @return 
#' an xts object with
#' Bid, Ask, Mid columns, 
#' or Open, Close, Adjusted columns, 
#' or Open, Close columns.
#' or Price column.
#' @author Lance Levenson, Brian Peterson, Garrett See
#' @note requires quantmod
#' @seealso 
#' \code{\link{buildSpread}}
#' \code{\link{synthetic.instrument}}
#' \code{\link{formatSpreadPrice}}
#' \code{\link{buildRatio}}
#' @examples
#' \dontrun{
#' currency("USD")
#' stock("SPY", "USD")
#' stock("DIA", "USD")
#' getSymbols(c("SPY","DIA"))
#' 
#' #can call with names of instrument/xts ojects
#' fSB <- fn_SpreadBuilder("SPY","DIA") 
#' fSB2 <- fn_SpreadBuilder(SPY,DIA) # or you can pass xts objects
#'
#' #assuming you first somehow calculated the ratio to be a constant 1.1
#' fSB3 <- fn_SpreadBuilder("SPY","DIA",1.1) 
#' head(fSB)
#'
#' # Call fn_SpreadBuilder with vector of 2 instrument names
#' # in 1 arg instead of using both prod1 and prod2.
#' fSB4 <- fn_SpreadBuilder(c("SPY","DIA"))
#' #download data and plot the closing values of a spread in one line
#' chartSeries(Cl(fn_SpreadBuilder(getSymbols(c("SPY","DIA")),auto.assign=FALSE)))
#' }
#' @import xts
#' @export
fn_SpreadBuilder <- function(prod1, prod2, ratio=1, currency='USD', from=NULL, 
    to=NULL, session_times=NULL, notional=TRUE,
    unique_method=c('make.index.unique','duplicated','least.liq','price.change'), 
    silent=FALSE, auto.assign=TRUE, env=.GlobalEnv, ...)
{
##TODO: allow for different methods for calculating Bid and Ask 
    dargs <- list(...)
    # dots may have args for getSymbols or for make_spread_id.
    # extract args that should be passed to make_spread_id
    msi.args <- list() #make_spread_id args
    for (arg in c("root", "format", "sep")) {
        if (!is.null(dargs[[arg]])) {
            msi.args[[arg]] <- dargs[[arg]]
            dargs[[arg]] <- NULL
        }
    }

    if (length(prod1) == 2 && missing(prod2)) {
        prod2 <- prod1[2]
        prod1 <- prod1[1]
    }

    unique_method<-unique_method[1]

    Data.1 <- NULL
    Data.2 <- NULL
    
    if (is.xts(prod1)) {
        Data.1 <- prod1
        prod1 <- deparse(substitute(prod1))
    }
    if (is.xts(prod2)) {
        Data.2 <- prod2
        prod2 <- deparse(substitute(prod2))
    }

    prod1.instr <- try(getInstrument(prod1, silent=TRUE))
    if (!is.instrument(prod1.instr) || inherits(prod1.instr,'try-error') 
        || !isTRUE(notional)) { 
        if (!silent && isTRUE(notional)) {
            warning(paste('could not find instrument ', prod1, 
                          '. Using multiplier of 1 and currency of ', 
                          currency, sep=''))
        }
        prod1.instr <- list(multiplier=1, currency=currency)
    }

    prod2.instr <- try(getInstrument(prod2, silent=TRUE))
    if (!is.instrument(prod2.instr) || inherits(prod2.instr, 'try-error') 
        || !isTRUE(notional)) {
        if (!silent && isTRUE(notional)) {
            warning(paste('could not find instrument ', prod2, 
                          '. Using multiplier of 1 and currency of ', 
                           currency, sep=''))
        }
        prod2.instr <- list(multiplier=1,currency=currency)
    }
    if (is.null(Data.1)) Data.1 <- try(get(as.character(prod1), pos=env),
                                           silent=TRUE) 
    if (is.null(Data.2)) Data.2 <- try(get(as.character(prod2), pos=env),
                                           silent=TRUE) 

    if (inherits(Data.1, "try-error") || (inherits(Data.2, "try-error"))) {
        gS.args <- list()
        if (!is.null(from)) gS.args$from <- from
        if (!is.null(to)) gS.args$to <- to
        gS.args$auto.assign=FALSE
        if (inherits(Data.1, 'try-error')) {
            Data.1 <- do.call("getSymbols", c(Symbols=prod1, gS.args, dargs))
        }
        if (inherits(Data.2, 'try-error')) {
            Data.2 <- do.call("getSymbols", c(Symbols=prod2, gS.args, dargs))
        }
    }

    
    if ( (all(has.Op(Data.1), has.Cl(Data.2)) && !(all(has.Op(Data.2), has.Cl(Data.2)))) || 
	(is.BBO(Data.1) && !is.BBO(Data.2)) ||
	(!(all(has.Op(Data.1), has.Cl(Data.2))) && (all(has.Op(Data.2), has.Cl(Data.2)))) ||
	(!is.BBO(Data.1) && is.BBO(Data.2)) ) {
        stop('prod1 and prod2 must be the same types of data (BBO,OHLC,etc.)')
    }
    
    if (is.null(from)) from <- max(index(first(Data.1)),index(first(Data.2)))
    if (is.null(to)) to <- min(index(last(Data.1)),index(last(Data.2))) 
    Data.1 <- Data.1[paste(from,to,sep="::")]
    Data.2 <- Data.2[paste(from,to,sep="::")]
    
    Mult.1 <- as.numeric(prod1.instr$multiplier) 
    Mult.2 <- as.numeric(prod2.instr$multiplier) 

    if (prod1.instr$currency != currency) {  
        #Cur.1 <- .get_rate(prod1.instr$currency,currency)
        Data.1 <- redenominate(prod1,currency)
    }
    if (prod2.instr$currency != currency) {  
        #Cur.2 <- .get_rate(prod2.instr$currency,currency)
        Data.2 <- redenominate(prod2,currency)
    }

    #Determine what type of data it is
    if (all(has.Op(Data.1), has.Cl(Data.1), has.Ad(Data.1))) {
      	M <- merge(Op(Data.1)[,1],Cl(Data.1)[,1],Ad(Data.1)[,1],Op(Data.2)[,1],
                   Cl(Data.2)[,1],Ad(Data.2)[,1])
	colnames(M) <- c("Open.Price.1","Close.Price.1","Adjusted.Price.1",
                     "Open.Price.2","Close.Price.2","Adjusted.Price.2")
    } else if(all(has.Op(Data.1), has.Cl(Data.1))) {
	M <- merge(Op(Data.1)[,1],Cl(Data.1)[,1],Op(Data.2)[,1],Cl(Data.2)[,1])
	colnames(M) <- c("Open.Price.1","Close.Price.1","Open.Price.2",
                     "Close.Price.2")
    } else if (is.BBO(Data.1)) {
	M <- merge(Data.1[,c( grep('Bid',colnames(Data.1),ignore.case=TRUE)[1], 
			grep('Ask',colnames(Data.1),ignore.case=TRUE)[1])],
		Data.2[,c(grep('Bid',colnames(Data.1),ignore.case=TRUE)[1],
			grep('Ask',colnames(Data.2),ignore.case=TRUE)[1])] )
      colnames(M) <- c("Bid.Price.1","Ask.Price.1","Bid.Price.2","Ask.Price.2")
    } else M <- merge(Data.1,Data.2)

    fn_split <- function(DF)
    {   
        DF.split <- split(DF,"days")
        ret <- NULL
        
        for(d in 1:length(DF.split))
        {
            tmp <- na.locf(DF.split[[d]])
            tmp <- na.omit(tmp)
            ret <- rbind(ret,tmp)   
        }
        #attr(attr(ret,"index"),"tzone") <- "GMT" # no longer needed?
        #attr(ret,".indexTZ") <- "GMT" # no longer needed?
	    colnames(ret) <- colnames(DF)
        ret
    }
    
    M <- fn_split(M)
	
    #can't subset times until after the merge
    if(!is.null(session_times)){
        #Data.1 <- Data.1[time.sub.GMT]
        #Data.2 <- Data.2[time.sub.GMT]
        M <- M[session_times]
    }
    
    if( all(has.Op(Data.1), has.Cl(Data.1)) ) {
      M$Open.Price.1 <- M$Open.Price.1 * Mult.1     # * Cur.1 
      M$Close.Price.1 <- M$Close.Price.1 * Mult.1   # * Cur.1
      M$Open.Price.2 <- M$Open.Price.2 * Mult.2     # * Cur.2
      M$Close.Price.2 <- M$Close.Price.2 * Mult.2   # * Cur.2
      
      open <- M$Open.Price.1 - M$Open.Price.2 * ratio
      close <- M$Close.Price.1 - M$Close.Price.2 * ratio

      Spread <- cbind(open,close)
      colnames(Spread) <- c('Open.Price','Close.Price')
      if (has.Ad(Data.1)) {
	    M$Adjusted.Price.1 <- M$Adjusted.Price.1 * Mult.1   # * Cur.1
	    M$Adjusted.Price.2 <- M$Adjusted.Price.2 * Mult.2   # * Cur.2
	    Spread$Adjusted.Price <- M$Adjusted.Price.1 - M$Adjusted.Price.2 * ratio
      }
      #Spread$Mid.Price <- (Spread$Open.Price + Spread$Close.Price) / 2
    } else
    if (is.BBO(Data.1) ) {
      M$Bid.Price.1 <- M$Bid.Price.1 * Mult.1 # * Cur.1 
      M$Ask.Price.1 <- M$Ask.Price.1 * Mult.1 # * Cur.1
      M$Bid.Price.2 <- M$Bid.Price.2 * Mult.2 # * Cur.2
      M$Ask.Price.2 <- M$Ask.Price.2 * Mult.2 # * Cur.2
      ##TODO: Expand this to work with multiple legs
      bid <- M$Bid.Price.1 - ratio * M$Ask.Price.2
      ask <- M$Ask.Price.1 - ratio * M$Bid.Price.2
      
      Spread <- cbind(bid,ask)
      names(Spread) <- c("Bid.Price","Ask.Price")
      Spread$Mid.Price <- (Spread$Bid.Price + Spread$Ask.Price) / 2
    } else {
    #univariate spread. 
      if (ncol(M) > 2) stop('Unrecognized column names.')
      Spread <- M[,1] - ratio * M[,2]
      colnames(Spread) <- paste(prod1,prod2,'Price',sep='.')
    }
##TODO: Test with symbols where each symbol has data on a day that the other one doesn't 
##TODO: Add a method that merges Data.1 and Data.2 with all=FALSE and use that index to subset
    switch(unique_method,
            make.index.unique = {Spread<-make.index.unique(Spread)},
            least.liq = {
                #determine the least liquid
                idx1 <- index(na.omit(getPrice(Data.1))) 
                idx2 <- index(na.omit(getPrice(Data.2)))
                if(length(idx1)<length(idx2)) idx<-idx1 else idx <- idx2
                
                #subset the Spread
                Spread <- Spread[idx]
            },
            duplicated = {
                Spread <- Spread[!duplicated(index(Spread))]  #this may still be useful for instrument with huge numders of observations 
            },
            price.change = {
                Spread <- Spread[which(diff(Spread$Mid.Price)!=0 | 
                                        diff(Spread$Bid.Price)!=0 | 
                                        diff(Spread$Ask.Price)!=0) ,]
                
            }
    )
    if (auto.assign) { #store the data in 
        msi.args$x <- c(prod1,prod2)
        id <- do.call("make_spread_id", msi.args) #can pass 'root' or 'format' through dots
        memberratio <- if(length(ratio) > 1) { 
            list(1,-as.numeric(ratio))
        } else c(1,-ratio)
        spread(id, currency=currency, members=msi.args$x, 
               memberratio=memberratio, defined.by='fn_SpreadBuilder')
        assign(id, Spread, pos=env)
        id
    } else Spread  
}


#' Construct a primary_id for a \code{spread} \code{instrument} from the 
#' primary_ids of its members
#'
#' @param x character vector of member primary_ids
#' @param root Optional character string of root_id to use.
#' @param format String indicating how to format the suffix_ids of the spread.  
#'    If \code{NULL} (the default), or \code{FALSE}, no formatting will be done.  
#'    See \code{\link{format_id}} for other accepted values for \code{format}
#' @param sep character string to separate root_id and suffix_id
#' @return character string that can be used as a primary_id for a 
#'    \code{\link{spread}} instrument
#' @author Garrett See
#' @seealso \code{\link{spread}}, \code{\link{build_spread_symbols}},  
#'    \code{\link{build_series_symbols}}
#' @examples
#' ids <- c('VX_aug1','VX_U11')
#' make_spread_id(ids, format='CY')
#' make_spread_id(ids, format=FALSE)
#' make_spread_id(c("VIX_JAN11","VIX_FEB11"),root='VX',format='CY')
#' @export
make_spread_id <- function(x, root=NULL, format=NULL, sep="_"){
#    if (length(x) != 2) stop("x must be a vector of length 2")
    if (is.null(root)) {
        root <- unlist(unique(sapply(x, parse_id)['root',]))
    }
    if (length(root) != 1) return(paste(format_id(x, format=format),collapse="."))
    suff <- paste(unlist(unique(sapply(x,parse_id)['suffix',])),collapse='.')
	#if (is.character(format)) suff <- paste(sapply(strsplit(suff,"\\.")[[1]], 
    #    format_id, format=format, parse='suffix'), collapse=".")
    if (!is.null(format) && is.character(format)) 
		suff <- paste(format_id(strsplit(suff,"\\.")[[1]], format=format, parse='suffix'), collapse=".")
    id <- paste(root,suff, sep=sep)
    return(make.names(id))
}


#' format the price of a synthetic instrument
#'
#' Divides the notional spread price by the spread multiplier and rounds prices 
#' to the nearest \code{tick_size}.
#' @param x xts price series
#' @param multiplier numeric multiplier (e.g. 1000 for crack spread to get 
#'   from $ to $/bbl)
#' @param tick_size minimum price change of the spread
#' @return price series of same length as \code{x}
#' @author Garrett See
#' @seealso
#' \code{\link{buildSpread}}, \code{\link{fn_SpreadBuilder}}
#' @export
formatSpreadPrice <- function(x,multiplier=1,tick_size=0.01) {
  x <- x / multiplier
  round( x / tick_size) * tick_size
}

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
# $Id: buildSpread.R 1638 2014-10-08 03:43:16Z gsee $
#
###############################################################################
