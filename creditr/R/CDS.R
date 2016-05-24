#' Build a \code{CDS} class object given the input about a CDS contract.
#' 
#' @name CDS
#'   
#' @param name is the name of the reference entity. Optional.
#' @param contract is the contract type, default SNAC
#' @param RED alphanumeric code assigned to the reference entity. Optional.
#' @param date is when the trade is executed, denoted as T. Default is
#'   \code{Sys.Date}. The date format should be in "YYYY-MM-DD".
#' @param spread CDS par spread in bps.
#' @param maturity date of the CDS contract.
#' @param tenor of contract. By default is set as 5
#' @param coupon quoted in bps. It specifies the payment amount from the
#'   protection buyer to the seller on a regular basis. The default is 100 bps.
#' @param recovery in decimal. Default is 0.4.
#' @param currency in which CDS is denominated.
#' @param notional is the amount of the underlying asset on which the payments
#'   are based. Default is 10000000, i.e. 10MM.
#'   
#' @return a \code{CDS} class object including the input informtion on the
#'   contract as well as the valuation results of the contract.
#'   
#' @examples
#' x <- CDS(date = as.Date("2014-05-07"), tenor = 5, spread = 50, coupon = 100) 

CDS <- function(name = NULL,
                contract = "SNAC",
                RED      = NULL,
                date     = Sys.Date(),
                spread   = NULL,
                maturity = NULL,
                tenor    = NULL,
                coupon   = 100,
                recovery = 0.4,
                currency = "USD",
                notional = 10000000){
  
  ## if all three of date, tenor and maturity are given as input,
  ## then we need to check if the three are compatible
  
  if(!(is.null(tenor) | is.null(maturity))){
    
    ## maturityshouldbe is the correct maturity given the input date and tenor
    
    maturityshouldbe <- add_dates(data.frame(date = date, tenor = tenor,
                                             currency = currency))$endDate
    
    ## check if the input maturity matches the correct should-be maturity
    
    stopifnot(maturity == maturityshouldbe) 
  }

  if ((is.null(spread))) stop("Please input spread")
 
  ## if maturity date is not given we use the tenor and vice-versa, to get dates
  ## using add_dates function. Results are stored in cdsdates
  
  if(is.null(maturity)){
    cdsDates <- add_conventions(add_dates(data.frame(date = as.Date(date),
                                                     tenor = tenor, currency = currency)))
  } else{
    if(is.null(tenor)){
      cdsDates <- add_conventions(add_dates(data.frame(date = as.Date(date),
                                                       maturity = as.Date(maturity),
                                                       currency = currency)))
    }

    ## if both are entered, we arbitrarily use one of them
    
    if((!is.null(tenor) & !is.null(maturity))){
      cdsDates <- add_conventions(add_dates(data.frame(date = as.Date(date),
                                                       maturity = as.Date(maturity),
                                                       currency = currency)))
    }
  }

  ## if entity name and/or RED code is not provided, we set it as NA
  
  if (is.null(name)) name <- "NA"
  if (is.null(RED)) RED <- "NA"

  ## if tenor is NULL, we determine the tenor using the maturity date
  
  if(is.null(tenor)){
    lt1 <- as.POSIXlt(as.Date(date, origin="1900-01-01"))
    lt2 <- as.POSIXlt(as.Date(maturity, origin="1900-01-01"))
    tenor <- as.numeric((lt2$year*12 + lt2$mon) - (lt1$year*12 + lt1$mon))/12
  }
  
  ## if maturity date is NULL, we set maturity date as the endDate, which
  ## obtained using add_dates.
  
  if(is.null(maturity)) maturity = as.Date(cdsDates$endDate)

  ## create object of class CDS using the data we extracted
  
  cds <- new("CDS",
             name     = name,
             contract = contract,
             RED      = RED,
             date     = date,
             maturity = maturity,
             tenor    = as.numeric(tenor),
             coupon   = coupon,
             recovery = recovery,
             currency = currency,
             notional = notional,
             spread   = spread)
  
  ## With a given spread is given, calculate principal and accrual

  df <- data.frame(date             = as.Date(cds@date),
                   spread           = spread,
                   coupon           = cds@coupon,
                   maturity         = as.Date(cds@maturity),
                   currency         = cds@currency,
                   recovery         = cds@recovery,
                   notional         = cds@notional,
                   stringsAsFactors = FALSE)
    
  cds@principal <- spread_to_upfront(x = df, notional = cds@notional, 
                                                  isPriceClean = TRUE)
    
  ## dirty upfront
    
  cds@upfront <- spread_to_upfront(x = df, notional = cds@notional, isPriceClean = FALSE)

  cds@accrual <- cds@upfront - cds@principal
  
  cds@spread.DV01 <- spread_DV01(df)
  cds@IR.DV01     <- IR_DV01(df) 
  cds@rec.risk.01 <- rec_risk_01(df)
  cds@pd          <- spread_to_pd(data.frame(spread = cds@spread,
                                             currency = cds@currency,
                                             recovery = cds@recovery,
                                             tenor = cds@tenor,
                                             date = cds@date))
  
  ## calculate the default exposure of a CDS contract based on the
  ## formula: Default Exposure: (1-Recovery Rate)*Notional - Principal
  
  cds@price <- (1 - cds@principal / notional) * 100
  
  ## return object with all the calculated data
  
  return(cds)
}