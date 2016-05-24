#' Calculate CS10
#' 
#' \code{CS10} calculates the change in upfront value when the spread rises by 
#' 10%, also known as the CS10 of a contract.
#' 
#' @param x data frame, contains all the relevant columns.
#' @param date.var character, column in x containing date variable.
#' @param currency.var character, column in x containing currency.
#' @param maturity.var character, column in x containing maturity date.
#' @param tenor.var character, column in x containing tenors.
#' @param spread.var character, column in x containing spread in basis points.
#' @param coupon.var character, column in x containing coupon rates in basis 
#'   points. It specifies the payment amount from the protection buyer to the 
#'   seller on an annual basis.
#' @param recovery.var character, column in x containing recovery rates. ISDA
#'   model standard recovery rate asscumption is 0.4.
#' @param notional.var character, column in x containing the amount of the 
#'   underlying asset on which the payments are based.
#' @param notional numeric, the notional amount for all pricing if there isn't a
#'   notional.var
#' @param recovery numeric, the recovery rate for all pricing if there isn't a 
#'   recovery.var
#'   
#' @return a vector containing the change in upfront in units of currency.var 
#'   when spread increase by 10%, for each corresponding CDS contract.
#'   
#' @examples 
#' x <- data.frame(date = as.Date(c("2014-04-22", "2014-04-22")),
#'                 currency = c("USD", "EUR"),
#'                 tenor = c(5, 5),
#'                 spread = c(120, 110),
#'                 coupon = c(100, 100),
#'                 recovery = c(0.4, 0.4),
#'                 notional = c(10000000, 10000000),
#'                 stringsAsFactors = FALSE)
#' CS10(x)

CS10 <- function(x,
                 date.var      = "date",
                 currency.var  = "currency",
                 maturity.var  = "maturity",
                 tenor.var     = "tenor",
                 spread.var    = "spread",
                 coupon.var    = "coupon",
                 recovery.var  = "recovery",
                 notional.var  = "notional",
                 notional      = 10000000,
                 recovery      = 0.4){
  
  ## check if certain variables are contained in x
  
  x <- check_inputs(x,
                    date.var      = date.var,
                    currency.var  = currency.var,
                    maturity.var  = maturity.var,
                    tenor.var     = tenor.var,
                    spread.var    = spread.var,
                    coupon.var    = coupon.var,
                    notional.var  = notional.var,
                    notional      = notional,
                    recovery.var  = recovery.var,
                    recovery      = recovery)
  
  x <- add_conventions(add_dates(x))
  
  CS10 <- rep(NA, nrow(x))
  
  for(i in 1:nrow(x)){
    
    ## extract currency specific interest rate data and date conventions using
    ## get_rates()
    
    rates.info <- get_rates(date = x$date[i], currency = x$currency[i])
    
    CS10[i] <- call_ISDA(x = x[i,], name = "CS10", rates.info = rates.info)
  }
  
  return(CS10)
}