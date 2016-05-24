#' Get rates from Markit
#' 
#' \code{download_markit} takes a data frame of dates and returns a data frame 
#' with the yields for different maturities. 
#' 
#' @param start date for gathering interest rates. Must be a Date type
#' @param end date for gathering interest rates. Must be a Date type
#' @param currency for which rates are being retrieved
#'   
#' @return data frame containing the rates from every day from start to end 
#'   dates. Note: the date in the output data frame does not refer to the rates 
#'   of that day but to the date on which the CDS is being priced. So the 
#'   corresponding rate is actually the rate of the previous day. Example: if 
#'   the column reads 2014-04-22, the corresponding rates are actually for 
#'   2014-04-21.
#'   
#' @import zoo
#' @import xts
#' @importFrom utils head tail
#'   
#' @seealso \code{\link{download_FRED}} \code{\link{build_rates}}
#'   
#' @examples
#' \dontrun{
#' download_markit(start = as.Date("2005-12-31"), end = as.Date("2006-01-04"), 
#'                 currency = "JPY")
#' }       

download_markit <- function(start, end, currency = "USD"){
  
  stopifnot(inherits(start, "Date"))
  stopifnot(inherits(end, "Date"))
  
  ## start date must be smaller than end date
  
  stopifnot(end >= start)
  
  ## Our loop will just append to x each time through.
  
  x <- NULL
  
  ## we let start date go back one day so that we get previous day's rate
  
  start <- start - 1
  
  while(end > start){
    
    ## use get_rates with the end date and keep reducing it by 1 day store this
    ## data in Rates
    
    Rates <- try(get_raw_markit(date = end, currency=currency))
    
    ## we use try so that if there is a date where rates are unavailable, it
    ## doesn't stop the function.
    
    if(is(Rates, "try-error")){
      end <- end - 1 
    } else {
      
      ## append all the data from the different dates where rates are available
      
      x <- rbind(x, data.frame(date = rep(end, length(as.numeric(as.character(Rates$rate)))),
                               currency = currency,
                               expiry = Rates$expiry,
                               rate = as.numeric(as.character(Rates$rate))))
      
      end <- end - 1      
    }    
  }
  
  ## convert factor type to character
  
  x$currency = as.character(x$currency)
  x$expiry   = as.character(x$expiry)
  
  ## convert data frame to xts object; use hte first column of the data frame as
  ## the time index of the xts object
  
  x.xts <- xts::xts(x[, -1], order.by = x[, 1])
  
  ## Here we use "xts" package because it is good at manipulating missing dates.
  ## We define here an empty zoo object for the merge() later.
  
  empty <- zoo::zoo(order.by = seq.Date(utils::head(zoo::index(x.xts), 1), 
                                        utils::tail(zoo::index(x.xts), 1), by = "days"))
  
  ## merge() is a cool command that can merge two zoo objects together. Whenever
  ## merge() finds a missing row of date in raw.data, it will use the previous
  ## day's rate to fill the blank of the empty row of date in the empty zoo
  ## object. This just satisfies our need, because according to ISDA Standard
  ## Model, we use previous business day's rate for holidays and weekends. And
  ## this is exactly carried out by merge() and na.locf(). For further
  ## explanation, see stackoverflow.
  
  data <- suppressWarnings(zoo::na.locf(merge(x.xts, empty)))
  
  ## convert xts object to data frame; use time index of xts object as the first
  ## date column of the data frame
  
  x <- data.frame(date = zoo::index(x.xts), zoo::coredata(x.xts))
  
  ## convert factor type to character or numeric type Notice that converting
  ## factor type to a numeric type is a little bit more tedious, cuz if you
  ## directly use as.numeric(x$rate) some information may be lost (see ?factor);
  ## so convert its factor level is more correct
  
  x$currency <- as.character(x$currency)
  x$expiry   <- as.character(x$expiry)
  x$rate     <- as.numeric(levels(x$rate))[x$rate]
  
  return(x)
}