#' Calculate IR.DV01
#' 
#' \code{IR_DV01} calculate the amount of change in upfront when there is a
#' 1/1e4 increase in interest rate for a data frame of CDS contracts.
#' 
#' @inheritParams CS10
#'   
#' @return a vector containing the change in upfront when there is a 1/1e4 
#'   increase in interest rate, for each corresponding CDS contract.
#'   
#' @examples 
#' x <- data.frame(date = c(as.Date("2014-04-22"), as.Date("2014-04-22")),
#'                 currency = c("USD", "EUR"),
#'                 tenor = c(5, 5),
#'                 spread = c(120, 110),
#'                 coupon = c(100, 100),
#'                 recovery = c(0.4, 0.4),
#'                 notional = c(10000000, 10000000),
#'                 stringsAsFactors = FALSE)
#' IR_DV01(x)

IR_DV01 <- function(x,
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
  
  IR.DV01 <- rep(NA, nrow(x))
  
  for(i in 1:nrow(x)){
    
    ## extract currency specific interest rate data and date conventions using
    ## get_rates()
    
    rates.info <- get_rates(date = x$date[i], currency = x$currency[i])
    
    IR.DV01[i] <- call_ISDA(x = x[i, ], name = "IR.DV01", rates.info = rates.info)
  }
  
  return(IR.DV01)
}