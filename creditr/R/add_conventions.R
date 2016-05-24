#' Return accounting conventions
#' 
#' \code{add_conventions} takes a data frame with a currency.var column and
#' returns the same data frame with eight other columns of accounting
#' conventions added to it.
#' 
#' @param x a data frame containing all necessary information
#' @param currency.var a character indicating the name of currency column
#'   
#' @return a data frame with eight more columns of accounting conventions:
#'   badDayConvention (a character indicating how non-business days are
#'   converted), mmDCC (the day count convention of the instruments),
#'   mmCalendars (any calendar adjustment for the CDS), fixedDCC (the day count 
#'   convention of the fixed leg), floatDCC (the day count convention of the
#'   floating leg), fixedFreq (the frequency of the fixed rate of swap being
#'   paid), floatFreq (the frequency of the floating rate of swap being paid)
#'   and swapCalendars (any calendar adjustment for swap rate)
#'   
#' @references 
#' \url{
#' http://www.cdsmodel.com/cdsmodel/assets/cds-model/docs/c-code%20Key%20Functions-v1.pdf
#' }
#' 
#' @examples 
#' x <- data.frame(date = c(as.Date("2014-05-06"), as.Date("2014-05-07")),
#'                 currency = c("USD", "JPY"))
#' add_conventions(x)

add_conventions <- function(x, currency.var = "currency"){

  stopifnot(currency.var %in% names(x))
  
  for(i in 1:nrow(x)){

    x$badDayConvention[i] <- "M"        ## how non-business days are converted
    x$mmDCC[i]            <- "ACT/360"  ## day count convention of the instruments
    x$floatDCC[i]         <- "ACT/360"  ## day count convention of the floating leg

    if(x[[currency.var]][i] == "USD"){
      
      x$mmCalendars[i]   <- "none"     ## calendar adjustment for the CDS
      x$fixedDCC[i]      <- "30/360"   ## day count convention of the fixed leg
      x$fixedFreq[i]     <- "6M"       ## frequency of the fixed rate of swap being paid
      x$floatFreq[i]     <- "3M"       ## frequency of the floating rate of swap being paid.
      x$swapCalendars[i] <- "none"     ## input for holiday files to adjust to business days.
    } else{
      
      if(x[[currency.var]][i] == "JPY"){
        
        x$mmCalendars[i]   <- "TYO"
        x$fixedDCC[i]      <- "ACT/365"
        x$fixedFreq[i]     <- "6M"
        x$floatFreq[i]     <- "6M"
        x$swapCalendars[i] <- "TYO"
      } else{
        
        if(x[[currency.var]][i] == "EUR"){
          
          x$mmCalendars[i]   <- "none"
          x$fixedDCC[i]      <- "30/360"
          x$fixedFreq[i]     <- "1Y"
          x$floatFreq[i]     <- "6M"
          x$swapCalendars[i] <- "none"
        }
      }
    }
  }
  return(x)
}