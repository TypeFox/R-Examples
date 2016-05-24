#' Get Rates from FRED
#' 
#' \code{download_FRED} returns the deposits and swap rates for the day input,
#' along with the date conventions for that specific currency. The source is
#' FRED.
#' 
#' @param start is the start date of the data frame we want
#' @param end is the end date of the data frame we want
#' @param currency is the three-letter currency code. As of now, it works for
#'   USD, EUR, and JPY. The default is JPY.
#'   
#' @return a data frame that contains the rates based on the ISDA pecifications
#'   
#' @import zoo
#' @importFrom utils head tail
#'   
#' @seealso \code{\link{download_markit}} \code{\link{build_rates}}
#'   
#' @examples
#' \dontrun{
#' download_FRED(start = as.Date("2003-12-31"), end = as.Date("2005-01-04"), 
#'               currency = "JPY")
#' }             

download_FRED <- function(start = as.Date("2004-01-01"), 
                          end = as.Date("2005-01-04"), currency = "JPY"){
  
  ## We make start date and end date go back one day, because we use the 
  ## previous business day's interest rate for trading date pricing. By making
  ## them go back one day and then calculate and then make the output go forward
  ## one day (so that it prints what we want), we save lots of trouble in the
  ## middle of the function of adjusting the day, because we want to adjust the
  ## day when we are calling the data from FRED
  
  start <- start - 1
  end <- end - 1
    
  ## define an empty data.frame for loop rbind later rate.complete.df is the
  ## return object
  
  rate.complete.df <- data.frame()
  
  ## time.intv is the time interval of the dates we want. it is a vector Date
  ## type
  
  time.intv <- seq(0, difftime(end, start))
  
  ## define expiry to be twelve months. FRED does not have data for expiry over
  ## a year
  
  expiry <- c("1M", "2M", "3M", "6M", "1Y")

  ## every time the loop run below, it goes to get a new expiry
  
  for(i in 1:5){
    
    ## Because the naming of FRED data for 1M~9M is different from that of
    ## 10M~12M, so we have to use a "if" here
    
    if(expiry[i] %in% c("1M", "2M", "3M", "6M")){
      FRED.symbol <- paste(currency, expiry[i], "TD156N", sep = "")
    } else{
      FRED.symbol <- paste(currency, "12MD156N", sep = "")
    } 
    
    ## raw.data is to get the raw data from FRED using getSymbols() from
    ## "quantmod" package. We have to use suppressWarnings here because
    ## getSymbols() often generates a series of warnings that are actually not
    ## warnings. See stackoverflow for furtner explanations. The
    ## suppressWarnings here is harmless. Be sure to set env = NULL because if
    ## you store the data in an environment, error will occur.
    
    raw.data <- suppressWarnings(quantmod::getSymbols(Symbols = FRED.symbol, 
                                            warnings = FALSE, src='FRED', 
                                            env = NULL))
    
    ## Here we use "xts" package because it is good at manipulating missing 
    ## dates. We define here an empty zoo object for the merge() later.
    
    empty <- zoo::zoo(order.by = seq.Date(utils::head(zoo::index(raw.data), 1), 
                                          utils::tail(zoo::index(raw.data), 1), by = "days"))
    
    ## merge() is a cool command that can merge two zoo objects together. the
    ## raw.data we get has missing rows, for example, 2004-01-03 is absent. So
    ## we define an empty row that has 2004-01-03 but no rate. we use merge() to
    ## merge raw.data and empty zoo object, then whenever merge() finds a
    ## missing row of date in raw.data, it will use the previous day's rate to
    ## fill the blank of the empty row of date in the empty zoo object. This
    ## just satisfies our need, because according to ISDA Standard Model, we use
    ## previous business day's rate for holidays and weekends. And this is
    ## exactly carried out by merge() and na.locf(). For further explanation,
    ## see stackoverflow. then we use [start+time.intv], where time.intv is
    ## defined earlier, to get the dates we want from the giant data frame
    
    data <- zoo::na.locf(merge(raw.data, empty))[start + time.intv]
    
    ## rate.partial.df is a data frame containing only one expiry type of rates.
    ## then later we use rbind to put it into rate.complete.df. here, we use
    ## date = index(data) to put xts object into a data frame with its date
    ## index an independent column.
    
    rate.partial.df <- data.frame(date = zoo::index(data), currency = currency,
                                  expiry = expiry[i], zoo::coredata(data), 
                                  row.names = NULL)
    
    ## A stupid thing about using getSymbols() is that when it reads from a .csv
    ## file, it always set header = TRUE, and we can't change it. so it always
    ## bring the nasty header of the csv file into our data frame. And to delete
    ## a header is tedious, as is shown below.
    
    names(rate.partial.df) <- NULL
    
    names(rate.partial.df) <- c("date", "currency", "expiry", "rate")
    
    ## then we rbind() the rate.partial.df and rate.complete.df
    
    rate.complete.df <- rbind(rate.complete.df, rate.partial.df)
  }
  
  ## since the rates.RData is in a descending date order, we sort 
  ## rate.complete.df here
  
  rate.complete.df <- rate.complete.df[order((rate.complete.df$date), 
                                             decreasing = TRUE),] 
  
  ## since FRED is in percentage while markit is in decimal, we devide the rates
  ## by 100
  
  rate.complete.df$rate <- rate.complete.df$rate / 100
  
  ## Last step, don't forget to add one day to all the dates in the data frame,
  ## because we use the previous business day's interest rate for the trading
  ## date pricing. If you don't understand this, see ?rates.
  
  rate.complete.df$date <- rate.complete.df$date + 1
  
  ## convert factor type to character type.
  
  rate.complete.df$currency <- as.character(rate.complete.df$currency)
  rate.complete.df$expiry   <- as.character(rate.complete.df$expiry)
  
  return(rate.complete.df)
}
