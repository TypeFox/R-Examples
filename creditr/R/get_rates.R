#' Get interest rates from rates.RData or the Markit website
#' 
#' \code{get_rates} returns the deposits and swap rates for the day input, along
#' with the date conventions for that specific currency. The day input should be
#' a weekday. If not, go to the most recent previous weekday.
#' 
#' @param date Trade date. The rates for a trade date T are published on T-1
#'   weekday. This date refers to the day on which we want the CDS to be priced,
#'   not the date for the interest rates as the interest rates will be used is
#'   the day before the trade date. Eg. If we are trying to find the rates used
#'   to price a CDS on 2014-04-22, it will return the rates of 2014-04-21
#' @param currency the three-letter currency code. As of now, it works for USD,
#'   EUR, and JPY. The default is USD.
#'   
#' @return a data frame that contains date (the CDS pricing date), 
#'   
#' @examples
#' \dontrun{
#' get_rates(as.Date("2014-05-07"), currency = "USD")
#' }

get_rates <- function(date, currency){
  
  stopifnot(inherits(date, "Date"))
  stopifnot(currency %in% c( "USD", "EUR", "JPY"))
  
  ## check rates.RData first, before get rates from the internet.
  
  x <- get("rates")[get("rates")$date == date & get("rates")$currency == currency, ]
  x$type <- as.character(substr(x$expiry, start = nchar(x$expiry), 
                                stop = nchar(x$expiry)))
  
  if(nrow(x) != 0){
    for(k in 1:length(x$type)){
      if(x$type[k] == "Y") x$type[k] <- "S"
    }
    
    x$type[which(x$expiry == "1Y")] <- "M"
  }
  
  ## if the wanted rates are not in rates.RData, then go get rates
  ## from the Markit website, using download_markit
  
  if(nrow(x) == 0){
    
    x <- download_markit(start = date, end = date, currency = currency)    
    x$type <- as.character(substr(x$expiry, start = nchar(x$expiry), 
                                  stop = nchar(x$expiry)))
    ## if maturity is 1Y, it is of type M
    
    if(nrow(x) != 0){
      for(k in 1:length(x$type)){
        if(x$type[k] == "Y") x$type[k] <- "S"
      }
      
      x$type[which(x$expiry == "1Y")] <- "M"
    }
    
    x$expiry <- as.character(x$expiry)
  }
  return(x)
}
