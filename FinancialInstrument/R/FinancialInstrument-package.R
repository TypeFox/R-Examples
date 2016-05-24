

#' @importFrom quantmod getQuote is.HLC is.OHLC OHLC is.BBO has.Bid has.Ask 
#'  has.Op has.Ad has.Trade has.Price getPrice Op Cl Ad has.Cl getSymbols   
#'  getSymbolLookup setSymbolLookup yahooQF has.Vo Vo getOptionChain
#' @importFrom TTR stockSymbols runSum
#' @importFrom zoo na.locf as.zoo coredata is.zoo index
#' @importFrom quantmod importDefaults
NULL

# for packages in Suggests:
globalVariables(c("timeSeries", "readHTMLTable", "%dopar%", "foreach"))

#' Construct, manage and store contract specifications for trading
#'
#' Transaction-oriented infrastructure for defining tradable instruments based 
#' on their contract specifications.  Construct and manage the definition of any 
#' asset class, including derivatives, exotics and currencies.  Potentially 
#' useful for portfolio accounting, backtesting, pre-trade pricing and other 
#' financial research.  Still in active development.
#'
#' The FinancialInstrument package provides a construct for defining and storing 
#' meta-data for tradable contracts (referred to as instruments, e.g., stocks, 
#' futures, options, etc.).  It can be used to create any asset class and 
#' derivatives, across multiple currencies.  
#'
#' FinancialInstrument was originally part of a companion package, blotter, that 
#' provides portfolio accounting functionality.  Blotter accumulates 
#' transactions into positions, then into portfolios and an account.  
#' FinancialInstrument is used to contain the meta-data about an instrument, 
#' which blotter uses to calculate the notional value of positions and the 
#' resulting P&L.  FinancialInstrument, however, has plenty of utility beyond 
#' portfolio accounting, and was carved out so that others might take advantage 
#' of its functionality.
#' 
#' As used here, 'instruments' are S3 objects of type 'instrument' or a subclass 
#' thereof that define contract specifications for price series for a tradable 
#' contract, such as corn futures or IBM common stock.  When defined as 
#' instruments, these objects are extended to include descriptive information 
#' and contract specifications that help identify and value the contract.
#' 
#' A simple example of an instrument is a common stock.  An instrument can be 
#' defined in brief terms with an identifier (e.g., "IBM").  Beyond the primary 
#' identifier, additional identifiers may be added as well and will work as 
#' 'aliases'.  Any identifier will do -- Bloomberg, Reuters-RIC, CUSIP, etc. -- 
#' as long as it's unique to the workspace. In addition, a stock price will be 
#' denominated in a currency (e.g., "USD") and will have a specific tick size 
#' which is the minimum amount that the price can be quoted and transacted in 
#' (e.g., $0.01).  We also define a 'multiplier' that is used when calculating 
#' the notional value of a position or transaction using a quantity and price 
#' (sometimes called a contract multiplier).  For a stock it's usually '1'.  
#' 
#' More care is needed when dealing with complex instruments, like futures.  
#' First, we have to define a future as a root contract.  This root is not 
#' tradable unto itself, but is used to generate a series of futures which are 
#' tradable and expire through time.  The root contract will provide an 
#' identifier (e.g., 'C' for the CME's corn contract), a denomination currency, 
#' a multiplier (one futures contract will cover multiple items) and a minimum 
#' tick size.  From that definition, a series of expiring contracts can be 
#' generated ("C_H08", "C_Z08", etc.) by specifying a suffix to be associated 
#' with the series, usually something like 'Z9' or 'Mar10' denoting expiration 
#' and year.  As you might expect, options are treated similarly.  The package 
#' also includes constructors for certain synthetic instruments, such as 
#' spreads.
#' 
#' FinancialInstrument doesn't try to exhaust the possibilities of attributes, 
#' so it instead allows for flexibility.  If you wanted to add an attribute to 
#' tag the exchange the instrument is listed on, just add it when defining the 
#' instrument (e.g., \code{future('CL', multiplier=1000, currency="USD", 
#' tick_size=.01, exchange="CME", description="Crude Light futures")}).  Or, as 
#' you can see, we've found it useful to add a field with more slightly more 
#' detail, such as \code{description='IBM Common Stock'}.  You can also add
#' attribute after the instrument has been created using 
#' \code{\link{instrument_attr}} as shown in the examples section below.
#'
#' Defining instruments can be tedious, so we've also included a CSV loader, 
#' \code{\link{load.instruments}}, in the package, as well as some functions
#' that will update instruments with data downloaded from the internet.
#' See, e.g., \code{\link{update_instruments.yahoo}}, 
#' \code{\link{update_instruments.TTR}}, 
#' \code{\link{update_instruments.morningstar}},
#' \code{\link{update_instruments.iShares}}. You can also update an instrument
#' using the details of another one with 
#' \code{\link{update_instruments.instrument}} which can be useful for creating
#' a new future_series from an expiring one.
#' 
#' Once you've defined all these instruments (we keep hundreds or thousands of 
#' them in our environments), you can save the instrument environment using
#' \code{\link{saveInstruments}}.  When you start a fresh R session, you 
#' can load your instrument definitions using \code{loadInstruments}.  We 
#' maintain an instrument.RData file that contains definitions for all 
#' instruments for which we have market data on disk.  
#' 
#' You may want to use \code{\link{setSymbolLookup.FI}} to define 
#' where and how your market data are stored so that 
#' \code{\link[quantmod]{getSymbols}} will work for you. 
#' 
#' FinancialInstrument's functions build and manipulate objects that are stored 
#' in an environment named ".instrument" at the top level of the package 
#' (i.e. "FinancialInstrument:::.instrument") rather than the global 
#' environment, \code{.GlobalEnv}.  Objects may be listed using 
#' \code{ls_instruments()} (or many other ls_* functions).
#'
#' We store instruments in their own environment for two reasons.  First, it 
#' keeps the user's workspace less cluttered and lowers the probability of 
#' clobbering something.  Second, it allows the user to save and re-use the 
#' \code{.instrument} environment in other workspaces.  Objects created with 
#' FinancialInstrument may be directly manipulated as any other object, but in 
#' our use so far we've found that it's relatively rare to do so.  Use the 
#' \code{\link{getInstrument}} function to query the contract specs of a 
#' particular instrument from the environment.
#' 
#' @name FinancialInstrument-package
#' @aliases FinancialInstrument-package FinancialInstrument
#' @docType package
#' @author
#' Peter Carl,
#' Brian G. Peterson,
#' Garrett See,
#'
#' Maintainer: G See \email{gsee000@@gmail.com}
#' @keywords package
#' @seealso
#' \code{\link[xts:xts-package]{xts}},
#' \code{\link[quantmod:quantmod-package]{quantmod}},
#' \href{https://r-forge.r-project.org/R/?group_id=316}{blotter},
#' \href{http://cran.r-project.org/web/packages/PerformanceAnalytics/index.html}{PerformanceAnalytics},
#' \href{https://r-forge.r-project.org/R/?group_id=1113}{qmao, and twsInstrument}
#' 
#' @examples
#' \dontrun{
#' # Construct instruments for several different asset classes
#' # Define a currency and some stocks
#' require("FinancialInstrument")
#' currency(c("USD", "EUR")) # define some currencies
#' stock(c("SPY", "LQD", "IBM", "GS"), currency="USD") # define some stocks
#' exchange_rate("EURUSD") # define an exchange rate
#' 
#' ls_stocks() #get the names of all the stocks
#' ls_instruments() # all instruments
#' 
#' getInstrument("IBM")
# update instruments with data from the web
#' update_instruments.yahoo(ls_stocks())
#' update_instruments.TTR(ls_stocks()) # doesn't update ETFs
#' update_instruments.masterDATA(ls_stocks()) # only updates ETFs
#' getInstrument("SPY")
#' 
#' ## Compare instruments with all.equal.instrument method
#' all.equal(getInstrument("USD"), getInstrument("USD"))
#' all.equal(getInstrument("USD"), getInstrument("EUR"))
#' all.equal(getInstrument("SPY"), getInstrument("LQD"))
#'
#' ## Search for the tickers of instruments that contain words
#' find.instrument("computer") #IBM
#' find.instrument("bond")  #LQD
#' 
#' ## Find only the ETFs; update_instruments.masterDATA added a "Fund.Type" field
#' ## to the ETFs, but not to the stocks
#' ls_instruments_by("Fund.Type") # all instruments that have a "Fund.Type" field
#' 
#' # build data.frames with instrument attributes
#' buildHierarchy(ls_stocks(), "Name", "type", "avg.volume")
#' 
#' ## before defining a derivative, must define the root (can define the underlying 
#' ## in the same step)
#' future("ES", "USD", multiplier=50, tick_size=0.25, 
#'     underlying_id=synthetic("SPX", "USD", src=list(src='yahoo', name='^GSPC')))
#' 
#' # above, in addition to defining the future root "ES", we defined an instrument 
#' # named "SPX".  Using the "src" argument causes setSymbolLookup to be called.
#' # Using the "src" arg as above is the same as 
#' # setSymbolLookup(SPX=list(src='yahoo', name='^GSPC'))
#' getSymbols("SPX") # this now works even though the Symbol used by 
#'                   # getSymbols.yahoo is "^GSPC", not "SPX"
#' 
#' ## Back to the futures; we can define a future_series
#' future_series("ES_U2", identifiers=list(other="ESU2"))
#' # identifiers are not necessary, but they allow for the instrument to be found 
#' # by more than one name
#' getInstrument("ESU2") #this will find the instrument even though the primary_id 
#'                       #is "ES_U2"
#' # can also add indentifiers later
#' add.identifier("ES_U2", inhouse="ES_U12")
#' 
#' # can add an arbitrary field with instrument_attr
#' instrument_attr("ES_U2", "description", "S&P 500 e-mini")
#' getInstrument("ES_U2")
#' 
#' option_series.yahoo("GS") # define a bunch of options on "GS"
#' # option root was automatically created
#' getInstrument(".GS")
#' # could also find ".GS" by looking for "GS", but specifiying type
#' getInstrument("GS", type='option')
#' 
#' # if you do not know what type of instrument you need to define, try
#' instrument.auto("ESM3")
#' getInstrument("ESM3")
#' instrument.auto("USDJPY")
#' getInstrument("USDJPY")
#' 
#' instrument.auto("QQQ") #doesn't work as well on ambigous tickers 
#' getInstrument("QQQ")
#'
#' # Some functions that make it easier to work with futures
#' M2C() # Month To Code
#' M2C()[5]
#' M2C("may")
#' C2M() # Code To Month
#' C2M("J")
#' C2M()[7]
#' MC2N("G") # Month Code to Numeric
#' MC2N("H,K,M")
#' 
#' parse_id("ES_U3")
#' parse_id("EURUSD")
#' 
#' next.future_id("ES_U2")
#' next.future_id("ZC_H2", "H,K,N,U,Z")
#' prev.future_id("CL_H2", 1:12)
#' 
#' sort_ids(ls_instruments()) # sort by expiration date, then alphabetically for 
#'                            # things that don't expire.
#' 
#' format_id("ES_U2", "CYY")
#' format_id("ES_U2", "CYY", sep="")
#' format_id("ES_U2", "MMMYY")
#' 
#' ## Saving the instrument environment to disk
#' tmpdir <- tempdir()
#' saveInstruments("MyInstruments.RData", dir=tmpdir)
#' rm_instruments(keep.currencies=FALSE)
#' ls_instruments() #NULL
#' loadInstruments("MyInstruments.RData", dir=tmpdir)
#' ls_instruments()
#' unlink(tmpdir, recursive=TRUE)
#' 
#' #build a spread:
#' fn_SpreadBuilder(getSymbols(c("IBM", "SPY"), src='yahoo'))
#' head(IBM.SPY)
#' getInstrument("IBM.SPY")
#' 
#' # alternatively, define a spread, then build it
#' spread(members=c("IBM", "GS", "SPY"), memberratio=c(1, -2, 1))
#' buildSpread("IBM.GS.SPY") #Since we hadn't yet downloaded "GS", buildSpread 
#'                           #downloaded it temporarily
#' chartSeries(IBM.GS.SPY)
#' 
#' ## fn_SpreadBuilder will return as many columns as it can 
#' ## (Bid, Ask, Mid, or Op, Cl, Ad), but only works on 2 instrument spreads
#' ## buildSpread works with any number of legs, but returns a single price column
#' 
#' getFX("EUR/USD", from=Sys.Date()-499) # download exchange rate from Oanda
#' 
#' IBM.EUR <- redenominate("IBM", "EUR") #price IBM in EUR instead of dollars
#' chartSeries(IBM, subset='last 500 days', TA=NULL)
#' addTA(Ad(IBM.EUR), on=1, col='red')
#'
#'}
NULL

# Copied from quantmod because it's not exported
has.Mid <- function (x, which = FALSE) {
    colAttr <- attr(x, "Mid")
    if (!is.null(colAttr)) 
        return(if (which) colAttr else TRUE)
    loc <- grep("Mid", colnames(x), ignore.case = TRUE)
    if (!identical(loc, integer(0))) 
        return(ifelse(which, loc, TRUE))
    ifelse(which, loc, FALSE)
}

# copied from quantmod:::convert.time.series since it is not exported
convert.time.series <- function (fr, return.class) {
    if ("quantmod.OHLC" %in% return.class) {
        class(fr) <- c("quantmod.OHLC", "zoo")
        return(fr)
    }
    else if ("xts" %in% return.class) {
        return(fr)
    }
    if ("zoo" %in% return.class) {
        return(as.zoo(fr))
    }
    else if ("ts" %in% return.class) {
        fr <- as.ts(fr)
        return(fr)
    }
    else if ("data.frame" %in% return.class) {
        fr <- as.data.frame(fr)
        return(fr)
    }
    else if ("matrix" %in% return.class) {
        fr <- as.data.frame(fr)
        return(fr)
    }
    else if ("its" %in% return.class) {
        if ("package:its" %in% search() || suppressMessages(require("its", 
            quietly = TRUE))) {
            fr.dates <- as.POSIXct(as.character(index(fr)))
            fr <- its::its(coredata(fr), fr.dates)
            return(fr)
        }
        else {
            warning(paste("'its' from package 'its' could not be loaded:", 
                " 'xts' class returned"))
        }
    }
    else if ("timeSeries" %in% return.class) {
        if ("package:timeSeries" %in% search() || suppressMessages(require("timeSeries", 
            quietly = TRUE))) {
            fr <- timeSeries(coredata(fr), charvec = as.character(index(fr)))
            return(fr)
        }
        else {
            warning(paste("'timeSeries' from package 'timeSeries' could not be loaded:", 
                " 'xts' class returned"))
        }
    }
}
