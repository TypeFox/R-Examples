################################################################################
##
## $Id: AllClasses.R 1317 2009-01-13 16:36:12Z enos $
##
## Class definitions for the portfolioSim package.
##
################################################################################

## Load hook for methods

.onLoad <- function(lib, pkg) require(methods)

## I want the below class to represent atomic values that might have
## attributes (like POSIXt, Date).  I _don't_ want things like lists
## to be acceptable here.  What's a more elegant way?

setOldClass(c("Date"))
setOldClass(c("POSIXt"))
setClassUnion("orderable", c("numeric","character","logical", "POSIXt", "Date"))

setClassUnion("listOrNull", c("list", "NULL"))

## Virtual classes (interfaces, loosely) for retrieving data.  Calling
## 'getData' on 'simPortfolioData' should yield an object of class
## 'simPortfolioInput'.

setClass("simSummaryInterface", "VIRTUAL")

setClassUnion("simSummaryInterfaceOrNull", c("simSummaryInterface",
                                             "NULL"))

setClass("simTradesInterface", "VIRTUAL")

setClass("simTrades",
         representation = representation(
           period     = "orderable",
           trades     = "trades")
         ,
         prototype = prototype(
           period     = numeric(0),
           trades     = new("trades")
           ),
         validity = function(object){

           if(!validObject(object@trades)){
             return("Invalid trades.")
           }
           else{
             return(TRUE)
           }
         }
         )

setClass("stiFromSignal",
         representation = representation(
           in.var   = "character",
           type     = "character",
           size     = "characterOrNumeric",
           sides    = "character",
           equity   = "numeric",
           target   = "environment",
           rebal.on = "orderable",
           trading.style = "character",
           chunk.usd = "numeric",
           turnover = "numeric"
           ),

         prototype = prototype(
           in.var   = "signal",
           type     = "equal",
           size     = "quintile",
           sides    = c("long", "short"),
           equity   = 50000000,
           rebal.on = numeric(0),
           trading.style = "immediate",
           chunk.usd = 50000,
           turnover = Inf
           ),
         
         contains = "simTradesInterface")

setClass("stiPresetTrades",
         representation = representation(
           periods    = "orderable",
           sim.trades = "list"
           ),
         prototype = prototype(
           periods    = "orderable",
           sim.trades = list()
           ),
         contains = "simTradesInterface",
         validity = function(object){

           if(!all(sapply(object@sim.trades, class) == "trades")){
             return("All objects in sim.trades slot must be class trades.")
           }

           if(!isTRUE(all.equal(length(object@periods),
                                length(object@sim.trades)))){
             return(paste("List of trades in sim.trades slot must have",
                          "same length as periods slot."))
           }
         }
         )


## Similarly, calling 'getData' on 'simMarketData' will return an object
## of class 'simMarketInput'.

setClassUnion("data.frameOrNull", c("data.frame","NULL"))

setClass("simDataInterface", "VIRTUAL")

setClass("simData",
         representation = representation(data = "data.frameOrNull"),
         
         validity = function(object){

           ## Data items currently required of the market data input
           ## interface are:
           ##
           ## period - non-NA period of the cross-section.
           ## 
           ## id - non-NA _unique_ identifier for each security.
           ##
           ## start.price - price at the start of the period.  May be NA.
           ##
           ## end.price - price at the end of the period.  May be NA.
           ##
           ## ret - total return for this period.  May be NA.

           required.cols <- c("period","id","start.price",
                              "end.price","ret", "volume",
                              "universe")
           
           if(!is.null(object@data)){
             if(!all(required.cols %in% names(object@data))){
               need.cols <- required.cols[! required.cols %in% names(object@data)]
               return(paste("The following columns must exist in data slot:",
                            paste(need.cols, collapse = ",")))
             }

             ## We enforce character id's on portfolio objects as well
             ## -- probably want to open up.  Only thing stopping me
             ## is I need to be careful about factors and merging
             ## (maybe factor is _not_ allowable, but pretty much
             ## anything else is?  Or factors are turned into
             ## characters?).
             
             if(!is.character(object@data$id))
               return("Data column 'id' must have class character")

             if(!is.numeric(object@data$start.price))
               return("Data column 'start.price' must have class numeric")

             if(!is.numeric(object@data$end.price))
               return("Data column 'end.price' must have class numeric")

             if(!is.numeric(object@data$ret))
               return("Data column 'ret' must have class numeric")

             if(!is.numeric(object@data$volume))
               return("Data column 'volume' must have class numeric")

             if(!is.logical(object@data$universe))
               return("Data column 'universe' must have class logical")


             ## Now check NA's and duplicates.

             if(any(is.na(object@data$id)))
               return("NA's not allowed in id column")

             if(any(duplicated(object@data$id)))
               return("values in id column must be unique")

             if(any(is.na(object@data$period)))
               return("NA's not allowed in period column")

             if(any(is.na(object@data$universe)))
               return("NA's not allowed in universe column")

             ## Finally, perform some checks on the data itself.

             if(any(object@data$start.price < 0, na.rm = TRUE) ||
                any(object@data$end.price < 0, na.rm = TRUE) ||
                any(object@data$volume < 0, na.rm = TRUE))
               return(paste("Negative values not allowed in start.price,",
                            "end.price, or volume columns"))

             if(any(object@data$ret < -1, na.rm = TRUE))
               return("Values for ret cannot be less than -1")
             
             return(TRUE)
           }
         }
         )

## There are two types of data associated with a period result.
##
## First, instant data captures a snapshot of the world at some point
## in time, but without any notion of duration.  We are concerned with
## instant data from the start and end instant of the result period.
## Examples of instant data are portfolios and exposures.
##
## Second, period data describes what occurs during the period.  All
## data that involves the passing of time, such as returns,
## contribution, and turnver, belongs in the period data category.
##
## The period slot of a single period result object is an orderable
## tag for the period and refers to an interval of time.  The
## instantData objects stored in the start.data and end.data slots
## include the specific start and end instants that define the period.

setClass("instantData",
         representation = representation(
           instant      = "orderable",
           equity.long  = "numeric",
           equity.short = "numeric",
           size.long    = "numeric",
           size.short   = "numeric",
           holdings     = "portfolio",
           exposure     = "exposure"),

         prototype = prototype(
           instant      = numeric(0),
           equity.long  = 0,
           equity.short = 0,
           size.long    = 0,
           size.short   = 0,
           holdings     = new("portfolio"),
           exposure     = new("exposure")          
           )
         )

setClass("periodData",
         representation = representation(
           period            = "orderable",
           turnover          = "numeric",
           universe.turnover = "numeric",
           performance       = "performance",
           contribution      = "contribution",
           trades            = "trades"),
         prototype = prototype(
           period            = numeric(0),
           turnover          = 0,
           universe.turnover = 0,
           performance       = new("performance"),
           contribution      = new("contribution"),
           trades            = new("trades"))
         )

setClass("simResultSinglePeriod",
         representation = representation(
           start.data   = "instantData",
           end.data     = "instantData",
           period.data  = "periodData"
           ),
         prototype = prototype(
           start.data   = new("instantData"),
           end.data     = new("instantData"),
           period.data  = new("periodData")
           )
         )

## I'm sure there's a better name than 'data' for the slot of most
## importance in 'simResult', but it will have to do for now.

## It's difficult to have all data around all the time, so this object
## might turn into a specification of where data is kept, as opposed
## to containing the data itself.

setClass("simResult",
         representation = representation(
           freq              = "numeric",
           data              = "list",
           errors            = "list",
           type              = "character",
           summary.interface = "simSummaryInterfaceOrNull"),
         prototype = prototype(
           freq              = 1,
           data              = list(),
           errors            = list(),
           type              = "basic",
           summary.interface = NULL
           )
         )

setClass("simResultCollection",
         representation = representation(
           data       = "list"
           ),
         prototype = prototype(
           data       = list()
           ),
         validity = function(object){
           if(!all(sapply(object@data, class) == "simResult")){
             return("All objects in data slots must have class 'simResult'")
           }
           if(length(object@data) > 0 && is.null(names(object@data))){
             return("List of 'simResult' objects in data slot must have names")
           }
           return(TRUE)
         }
         )

## Some simple data connector extensions.

## This portfolio data interface gets signal data from a data frame
## and creates new portfolios from that signal data.

setClass("sdiDf",
         representation = representation(
           data       = "data.frame"),
         prototype = prototype(
           data       = data.frame()),
         contains = "simDataInterface")


setClass("portfolioSim",
         representation = representation(
           periods          = "data.frame",
           freq             = "numeric",

           ## Interface objects

           data.interface   = "simDataInterface",
           trades.interface = "simTradesInterface",
           summary.interface = "simSummaryInterfaceOrNull",

           ## Starting position for holdings
           
           start.holdings   = "portfolio",

           ## Fill parameters

           fill.volume.pct  = "numeric",
           
           ## Where will we specify exp.vars and things like this?
           ## For now, placing at the top level.

           exp.var          = "character",
           contrib.var      = "character",

           out.loc          = "character",
           out.type         = "character"),
         
         prototype = prototype(
           periods          = data.frame(period = numeric(0),
                                         start  = numeric(0),
                                         end    = numeric(0)),
           freq             = 1,

           data.interface   = new("sdiDf"),
           trades.interface = new("stiFromSignal"),
           
           start.holdings   = new("portfolio"),
           fill.volume.pct  = 15,
           exp.var          = character(0),
           contrib.var      = character(0),
           out.loc          = character(0),
           out.type         = "default"
           ),
         
         validity = function(object){
           period.cols <- c("period","start","end")
           
           if(!all(period.cols %in% names(object@periods))){
             return(paste("Columns required in periods data frame are",
                          paste(period.cols, collapse = " ")))
           }
           
           if(!is(object@periods$period, "orderable")){
             return("Period column of data frame in periods slot is not orderable")
           }

           if(!is(object@periods$start, "orderable")){
             return("Start column of data frame in periods slot is not orderable")
           }

           if(!is(object@periods$end, "orderable")){
             return("End column of data frame in periods slot is not orderable")
           }

           if(length(object@freq) != 1){
             return("freq must be a numeric vector of length 1")
           }

           valid.types <- c("basic", "detail", "contributions", "exposures",
                            "trades", "portfolio", "default", "all", "lean")

           if(!all(object@out.type %in% valid.types)){
             return("Invalid out.type found.")
           }
            
           if(any(duplicated(object@periods))){
             return("Each period should appear only once.")
           }

           return(TRUE)
         }
         )
