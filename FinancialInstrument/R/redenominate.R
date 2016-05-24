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
# $Id: redenominate.R 1638 2014-10-08 03:43:16Z gsee $
#
###############################################################################

#' get an exchange rate series
#'
#' Try to find exchange rate data in an environment, inverting if necessary.
#' 
#' @param ccy1 chr name of 1st currency
#' @param ccy2 chr name of 2nd currency
#' @param env environment in which to look for data.
#' @return xts object with as many columns as practicable.
#' @author Garrett See
#' @seealso 
#' \code{\link{buildRatio}}
#' \code{\link{redenominate}}
#' @examples
#'
#' \dontrun{
#' EURUSD <- getSymbols("EURUSD=x",src='yahoo',auto.assign=FALSE)
#' USDEUR <- .get_rate("USD","EUR")
#' head(USDEUR)
#' head(EURUSD)
#' }
#' @rdname get_rate
.get_rate <- function(ccy1, ccy2, env=.GlobalEnv) {
    rsym <- NA
    invert <- FALSE
    if (exists(paste(ccy1, ccy2, sep=""),where=env)) {
        rsym <- paste(ccy1, ccy2, sep="")
    } else if (exists(paste(ccy2, ccy1, sep=""),where=env)) {
        rsym <- paste(ccy2, ccy1, sep="")
        invert = TRUE
    } else if (exists(paste(ccy1, ccy2, sep="."),where=env)) {
        rsym <- paste(ccy1, ccy2, sep=".")
    } else if (exists(paste(ccy2, ccy1, sep="."),where=env)) {
        rsym <- paste(ccy2, ccy1, sep=".")
        invert = TRUE
    } else if (exists(paste(ccy1, ccy2, sep="/"),where=env)) {
        rsym <- paste(ccy1, ccy2, sep="/")
    } else if (exists(paste(ccy2, ccy1, sep="/"),where=env)) {
        rsym <- paste(ccy2, ccy1, sep="/")
        invert = TRUE
    }
    rate <- try(get(rsym,pos=env),silent=TRUE)
    if (inherits(rate,'try-error')) 
        stop(paste('Could not find exchange rate for ', ccy1, 
                    ' and ', ccy2, ' in ', deparse(substitute(env)), sep=''))
    rsym <- paste(substr(rsym,1,3), substr(rsym,nchar(rsym)-2,nchar(rsym)),sep="")
    if (invert) {
        rate <- 1/rate  #inverting will reverse High/Low and Bid/Ask
        rsym.inv <- paste(substr(rsym,4,6),substr(rsym,1,3),sep="")
        if (is.OHLC(rate)) {
            cn <- colnames(rate)
            hc <- grep('High',cn,ignore.case=TRUE)
            lc <- grep('Low',cn,ignore.case=TRUE)
            cn[hc] <- gsub('High','Low',cn[hc])
            cn[lc] <- gsub('Low','High',cn[lc])
            colnames(rate) <- gsub(rsym,rsym.inv,cn)
            rate <- OHLC(rate)
        } else if (is.BBO(rate)) {
            cn <- colnames(rate)
            bc <- grep('Bid',cn,ignore.case=TRUE)
            ac <- grep('Ask',cn,ignore.case=TRUE)
            cn[bc] <- gsub('Bid','Ask',cn[bc])
            cn[ac] <- gsub('Ask','Bid',cn[ac])
            colnames(rate) <- gsub(rsym,rsym.inv,cn)
            tmprate <- rate[, c(has.Bid(rate,1),has.Ask(rate,1))]
            rate <- cbind(tmprate,rate[,-c(has.Bid(rate,1),has.Ask(rate,1))])
        }
    }
    rate
}

#' Extract a single row from each day in an xts object
#' @param x xts object of sub-daily data.
#' @param EOD_time time of day to use.
#' @return xts object with daily scale.
#' @author Garrett See
#' @seealso quantmod:::to.daily, quantmod:::to.period
#' @export
#' @rdname to_daily
.to_daily <- function(x, EOD_time="15:00:00") {
    x <- do.call(rbind, lapply(split(x[paste("T00:00:00/T",EOD_time,sep="")],'days'),'last'))
    xts(x, order.by=as.Date(paste(index(x))))
}

#' construct price ratios of 2 instruments 
#'
#' Calculates time series of ratio of 2 instruments using available data. 
#' Returned object will be ratios calculated using Bids, Asks, and Mids, or Opens, Closes, and Adjusteds.
#'
#' \code{x} should be a vector of 2 instrument names. An attempt will be made to \code{get} the data
#' for both instruments.  If there are no xts data stored under either of the names, it will try to 
#' return prebuilt data with a call to \code{\link{.get_rate}}.
#'
#' If the data are not of the same frequency, or are not of the same type (OHLC, BBO, etc.)
#' An attempt will be made to make them compatible.  Preference is given to the first leg.
#'
#' If the data in \code{x[1]} is daily or slower and the data in \code{x[2]} is intraday
#  then the intraday data in \code{x[2]} will become univariate
#' (e.g. if you give it daily OHLC and intraday Bid Ask Mid, it will use all of 
#' the OHLC columns of \code{x[1]} and only the the End of Day Mid price of the BAM object.
#'
#' If the data in \code{x[1]} is intraday, and the data in \code{x[2]} is daily or slower,
#' for each day, the previous closing value of \code{x[2]} will be filled forward with \code{na.locf}
#'
#' @param x vector of instrument names. e.g. c("SPY","DIA")
#' @param env environment where xts data is stored
#' @param silent silence warnings?
#' @return 
#' An xts object with columns of
#' Bid, Ask, Mid
#' OR
#' Open, Close, Adjusted
#' OR
#' Price
#' @author Garrett See
#' @seealso
#' \code{\link{redenominate}}
#' \code{\link{buildSpread}}
#' \code{\link{fn_SpreadBuilder}}
#' @examples
#'
#' \dontrun{
#' syms <- c("SPY","DIA")
#' getSymbols(syms)
#' rat <- buildRatio(syms)
#' summary(rat)
#' }
#' @export
buildRatio <- function(x,env=.GlobalEnv, silent=FALSE) {
    if (length(x) != 2 || !is.character(x)) {
        stop('Please provide vector of 2 instruments names')
    }
    rat.sym <- paste(x[1],x[2],sep="")
    x1 <- try(get(x[1],pos=env),silent=TRUE)
    x2 <- try(get(x[2],pos=env),silent=TRUE)
    if (inherits(x1,'try-error') || inherits(x2,'try-error')) {
        #maybe we can get the ratio directly
        if (!silent) warning(paste('Nothing to build. Returning data found in', deparse(substitute(env)),'if any.'))
        return(.get_rate(x[1],x[2],env))  
    }
    #!#---#!# 
    Bi <- #This, or Bid, should be exported from quantmod
    function (x) 
    {
        if (has.Bid(x)) 
            return(x[,has.Bid(x,1)])            
            #return(x[, grep("Bid", colnames(x), ignore.case = TRUE)])
        stop("subscript out of bounds: no column name containing \"Bid\"")
    }
    As <- #This, or Ask, should be exported from quantmod
    function (x) 
    {
        if (has.Ask(x)) 
            return(x[,has.Ask(x,1)])    
            #return(x[, grep("Ask", colnames(x), ignore.case = TRUE)])
        stop("subscript out of bounds: no column name containing \"Ask\"")
    }
    
    Mid <- #This should be exported from quantmod
    function (x) 
    {
        if (has.Mid(x)) 
            return(x[,has.Mid(x,1)])            
            #return(x[, grep("Mid", colnames(x), ignore.case = TRUE)])
        stop("subscript out of bounds: no column name containing \"Mid\"")
    }
    #!#---#!#

    instr2 <- NULL
    instr1 <- try(getInstrument(x[1],silent=TRUE))
    if (is.instrument(instr1)) instr2 <- try(getInstrument(x[2],silent=TRUE))
    if (is.instrument(instr2)) {
        mult1 <- as.numeric(instr1$multiplier)
        mult2 <- as.numeric(instr2$multiplier)
    } else mult1 <- mult2 <- 1
    mrat <- mult1 / mult2

    if (is.OHLC(x1) && is.OHLC(x2)) {
        op <- Op(x1)[,1] / Op(x2)[,1] * mrat
        cl <- Cl(x1)[,1] / Cl(x2)[,1] * mrat
        if (!has.Ad(x1)) x1$Adjusted <- Cl(x1)[,1]
        if (!has.Ad(x2)) x2$Adjusted <- Cl(x2)[,1]
        ad <- Ad(x1)[,1] / Ad(x2)[,1] * mrat
        rat <- cbind(op,cl,ad)
        colnames(rat) <- paste(rat.sym, c("Open","Close","Adjusted"),sep='.')
    } else if (is.BBO(x1) && is.BBO(x2)) {
        bid <- Bi(x1)[,1]/As(x2)[,1] * mrat
        ask <- As(x1)[,1]/Bi(x2)[,1] * mrat
        if (has.Mid(x1) && has.Mid(x2)) {
            mid <- Mid(x1)[,1] / Mid(x2)[,1] * mrat 
        } else {
            mid <- ((Bi(x1)[,1]+As(x1)[,1])/2) / ((Bi(x2)[,1]+As(x2)[,1])/2) * mrat
        }
        rat <- cbind(bid,ask,mid)
        colnames(rat) <- paste(rat.sym,c('Bid','Ask','Mid'),sep='.')
    } else if (NCOL(x1) == 1 && NCOL(x2) == 1) {
        rat <- x1 / x2 * mrat #coredata(x1) / coredata(x2)
    } else if (periodicity(x1)$frequency >= 86400) { 
        #if daily or slower use OHLC and Mid
        if (is.OHLC(x1)) { #If first leg is.OHLC, 2nd leg will be univariate
            div <- if (NCOL(x2) == 1) {
                        x2
                   } else if (has.Mid(x2)) {
                        Mid(x2)[,1]
                   } else getPrice(x2)
            rat <- mrat * x1[,1] / div
            if (NCOL(x1) > 1) {
                for (i in 2:NCOL(x1)) {
                    rat <- cbind(rat, mrat * x1[,i]/div)
                }
            }
        } else if (is.OHLC(x2)) { #1st leg will be univariate
            num <- if (NCOL(x1) == 1){
                        x1
                    } else if (has.Mid(x1)) {
                        Mid(x1)[,1]
                    } else getPrice(x1)
            rat <- mrat * num / x2[,1]
            if (NCOL(x2) > 1) {
                for (i in 2:NCOL(x2)) {
                    rat <- cbind(rat, mrat * num/x2[,i])
                }
                colnames(rat) <- colnames(x2)                  
            }  
        }    
    } else if (periodicity(x1)$frequency < 86400) {
        #if intraday, use BAM and Cl
        if (is.BBO(x1)) { #1st leg is.BBO, 2nd leg will be univariate
            div <- if (NCOL(x2) == 1) {
                        x2
                    } else if (has.Cl(x2)) {
                        Cl(x2)[,1]
                    } else if (has.Ad(x2)) {
                        Ad(x2)[,1]    
                    } else getPrice(x2)[,1]
            rat <- mrat * x1[,1] / div
            if (NCOL(x1) > 1) {            
                for (i in 2:NCOL(x1)) {
                    rat <- cbind(rat, mrat * x1[,i]/div)
                }
            }
        } else if (is.BBO(x2)) { #1st leg will be univariate
            num <- if (NCOL(x1) == 1) {
                        x1
                    } else if (has.Cl(x1)) {
                        Cl(x1)[,1]
                    } else if (has.Ad(x1)) {
                        Ad(x1)[,1]
                    } else getPrice(x1)[,1]
            rat <- mrat * num / x2[,1]
            if (NCOL(x2) > 1){
                for (i in 2:NCOL(x2)) {
                    rat <- cbind(rat, mrat * num/x2[,i]) 
                }
            }
        }

    } else stop("I'm not programmed to handle this yet.")

    if (NCOL(rat) == 1)
        colnames(rat) <- paste(rat.sym,'price',sep='.')
    rat
}

#' Redenominate (change the base of) an instrument
#'
#' Redenominate (change the base of) an instrument
#'
#' If \code{old_base} is not provided, \code{x} must be the name of an 
#' instrument (or an object with the name of a defined instrument) so that the
#' currency attribute of the instrument can be used.  Otherwise, \code{old_base}
#' must be provided.
#'
#' If you want to convert to JPY something that is denominated in EUR,
#' you must have data for the EURJPY (or JPYEUR) exchange rate. If you don't have
#' data for EURJPY, but you do have data for EURUSD and USDJPY, 
#' you could \code{redenominate} to USD, then \code{redenominate} to EUR, 
#' but this function is not yet smart enough to do that for you.
#'
#' See the help for buildRatio also. 
#'
#' @param x can be either an xts object or the name of an instrument.
#' @param new_base change the denomination to this; usually a currency.
#' @param old_base what is the current denomination? 
#' @param EOD_time If data need to be converted to daily, this is the time of day to take the observation.
#' @param env environment that holds the data
#' @param silent silence warnings?
#' @return xts object, with as many columns as practicable, that represents the value of an instrument in a different currency (base).
#' @author Garrett See
#' @note this does not yet define any instruments or assign anything.
#' @seealso 
#' \code{\link{buildRatio}}
#' @examples
#'
#' \dontrun{
#' require(quantmod)
#' EURUSD <- getSymbols("EURUSD=x",src='yahoo',auto.assign=FALSE)
#' GLD <- getSymbols("GLD", src='yahoo', auto.assign=FALSE)
#' GLD.EUR <- redenominate(GLD,"EUR","USD") #can call with xts object
#'
#' currency("USD")
#' stock("GLD","USD")
#' GLD.EUR <- redenominate('GLD','EUR') #can also call with instrument name
#' }
#' @export
redenominate <- function(x, new_base='USD', old_base=NULL, EOD_time='15:00:00', env=.GlobalEnv, silent=FALSE) {
#TODO: create an instrument with currency=new_base.
    if (is.xts(x)) {
         Symbol <- deparse(substitute(x))
    } else Symbol <- x
    if (is.character(Symbol)) {
        instr <- try(getInstrument(Symbol,silent=TRUE))
        if (!is.instrument(instr)) {
            if (is.null(old_base)) stop(paste("If old_base is not provided, ", Symbol, ' must be defined.', sep=""))
            mult <- 1        
        } else {
            if (is.null(old_base)) old_base <- instr$currency
            mult <- as.numeric(instr$multiplier)    
        }
        if (is.character(x)) x <- get(Symbol,pos=env)
    }
    idxx <- index(x)
    #Now figure out the exchange rate
    #First assume that both bases are currencies, and look for an exchange rate
    if (!identical(new_base, old_base)) {
        rate <- try(.get_rate(new_base,old_base,env),silent=TRUE) #try with formats like EURUSD, EUR.USD, EUR/USD, and their inverses
        if (inherits(rate,'try-error')) {
            rate <- buildRatio(c(old_base, new_base), env=env) #maybe it's not FX
        }
    } else rate <- xts(rep(1L, nrow(x)), index(x))
    
    #!#---#!# Define function we'll need
    #This should be exported from quantmod
    Mid <- function (x) {
        if (has.Mid(x)) 
            return(x[,has.Mid(x,1)])            
            #return(x[, grep("Mid", colnames(x), ignore.case = TRUE)])
        stop("subscript out of bounds: no column name containing \"Mid\"")
    }
    #!#---#!#

    #Now we have data in x that needs to be multilied by data in rate.
    #First make sure they are the same periodicity

    #If you have daily data for x and intraday data for rate
    #convert rate to periodicity of x
    if (periodicity(x)$frequency >= 86400 && periodicity(rate)$frequency < 86400) { #x frequency is daily or lower, but rate freq is intraday
        if (is.OHLC(rate) || NCOL(rate) == 1) {
            rate <- to.period(rate, periodicity(x)$units)
        } else if(is.BBO(rate)) {
            if (periodicity(x)$scale == 'daily') {
                rate <- .to_daily(rate, EOD_time) #This doesn't make OHLC, the rest do.
            } else rate <- to.period(Mid(rate)[,1], periodicity(x)$units) 
        } else rate <- to.period(getPrice(rate)[,1], periodicity(x)$units)
    }

    # If you have intraday data for x and daily data for rate
    # use the daily rate for all rows of each day.
    if (periodicity(x)$frequency < 86400 && periodicity(rate)$frequency >= 86400) {
        df <- cbind(x, rate, all=TRUE)
        df <- na.locf(df,na.rm=TRUE)
        x <- df[, 1:NCOL(x)]
        rate <- df[, (NCOL(x)+1):NCOL(df)]
    }

    ff <- merge(rate,x,all=FALSE)
    ff <- na.omit(ff)
    ff <- ff[idxx]
    rate <- ff[,1:NCOL(rate)]
    x <- ff[,(NCOL(rate)+1):NCOL(ff)]

    tmpenv <- new.env()
    rsym <- new_base
    assign(rsym,rate,pos=tmpenv)
    assign(Symbol,x,pos=tmpenv)
    
    buildRatio(c(Symbol,rsym),env=tmpenv, silent=TRUE) / mult
#TODO: colnames
#TODO: auto.assign
}


#dailyConvertFX <- function(xts_obj, rate, prefer=NULL, EOD_time="11:00:00", verbose=TRUE) {
#    #to convert a EUR denominated asset from EUR to USD, rate=EURUSD
#    #DAX closes at 11:45 EDT or 10:45 Chicago time    
#    #FRED data is noon EDT or 11:00:00 Chicago time.
#    if (periodicity(xts_obj)$scale != "daily") stop('xts_obj must be daily')
#    rate <- getPrice(rate, prefer=prefer)    
#    tmpdt <- as.Date(index(rate[1:2,]))
#    if (tmpdt[1] == tmpdt[2]) { #intraday data
#        if (verbose) warning('converting rate to daily')
#        rate <- .to_daily(rate, EOD_time)        
#        rate <- rate[paste(start(xts_obj), end(xts_obj), sep="/")]
#    }   
#    df <- cbind(rate, xts_obj, all=TRUE)
#    df <- df[paste(max(start(rate),start(xts_obj)), '::', sep="")]
#    if (verbose && (NROW(df) < NROW(xts_obj))) warning('Data removed where rate was missing')
#    as.vector(df[,1]) * df[,2:(NCOL(df))]
#}


