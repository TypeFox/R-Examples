#' Calculate Spread with a Given Upfront
#' 
#' \code{upfront_to_spread} calculates conventional spread using the upfront 
#' or ptsUpfront values.
#' 
#' @param x data frame, contains all the relevant columns.
#' @param date.var character, column in x containing date variable.
#' @param currency.var character, column in x containing currency.
#' @param maturity.var character, column in x containing maturity date.
#' @param tenor.var character, column in x containing tenors.
#' @param coupon.var character, column in x containing coupon rates in basis 
#'   points. It specifies the payment amount from the protection buyer to the 
#'   seller on an annual basis.
#' @param recovery.var character, column in x containing recovery rates. ISDA
#'   model standard recovery rate asscumption is 0.4.
#' @param recovery numeric, the recovery rate for all pricing if there isn't a 
#'   recovery.var
#' @param points.var character name of points Upfront column
#' @param isPriceClean a boolean variable indicating whether the upfront is clean or dirty
#' @param notional numeric variable indicating the notional value of the CDS contract
#' @param payAccruedAtStart whether pay at start date the accrual amount
#' @param upfront.var is the character name of upfront column
#' @param payAccruedOnDefault whether pay in default scenario the accrual amount
#' 
#' @return a numeric indicating the spread.

upfront_to_spread <- function(x, 
                   currency.var = "currency", 
                   date.var = "date",
                   coupon.var = "coupon",
                   tenor.var = "tenor",
                   maturity.var = "maturity",
                   recovery.var = "recovery",
                   upfront.var = "upfront",
                   points.var = "ptsUpfront",
                   
                   isPriceClean = FALSE,
                   notional = 1e7,
                   payAccruedAtStart = FALSE,
                   payAccruedOnDefault = TRUE){
  
  x[[currency.var]] <- as.character(x[[currency.var]])
  
  if (is.null(x[[upfront.var]]) & is.null(x[[points.var]]))
    stop("Please input upfront or pts upfront")
  
  ## You must provide either a maturity or a tenor, but not both.
  
  stopifnot(!(is.null(x[[maturity.var]]) & is.null(x[[tenor.var]])))
  stopifnot(is.null(x[[maturity.var]]) | is.null(x[[tenor.var]]))

  
  if (is.null(x[[points.var]])) {
    x[[points.var]] <- x[[upfront.var]] / notional
  } else {
    payAccruedAtStart <- TRUE
  }
  
  x <- add_conventions(add_dates(x))
  
  spread <- rep(NA, nrow(x))
  
  for(i in 1:nrow(x)){
    
    rates.info <- get_rates(date = as.Date(x[[date.var]][i]),
                           currency = as.character(x[[currency.var]][i]))
    
    spread[i] <- .Call('calcCdsoneSpread',
                         baseDate_input = separate_YMD(x$baseDate[i]),
                         types = paste(as.character(rates.info$type), collapse = ""),
                         rates = as.numeric(as.character(rates.info$rate)),
                         expiries = as.character(rates.info$expiry),
                         mmDCC = x$mmDCC[i],
                         
                         fixedSwapFreq = x$fixedFreq[i],
                         floatSwapFreq = x$floatFreq[i],
                         fixedSwapDCC = x$fixedDCC[i],
                         floatSwapDCC = x$floatDCC[i],
                         badDayConvZC = x$badDayConvention[i],
                         holidays = x$swapCalendars[i],
                         
                         todayDate_input = separate_YMD(x[[date.var]][i]),
                         valueDate_input = separate_YMD(x$valueDate[i]),
                         benchmarkDate_input = separate_YMD(x$startDate[i]),
                         startDate_input = separate_YMD(x$startDate[i]),
                         endDate_input = separate_YMD(x$endDate[i]),
                         stepinDate_input = separate_YMD(x$stepinDate[i]),
                         
                         couponRate_input = x[[coupon.var]][i] / 1e4,
                         payAccruedOnDefault_input = payAccruedOnDefault,
                         
                         dccCDS = "ACT/360",
                         dateInterval = "Q",
                         stubType = "F",
                         badDayConv_input = "F",
                         calendar_input = "None",
                         
                         upfrontCharge_input = x[[points.var]][i],
                         recoveryRate_input = x[[recovery.var]][i],
                         payAccruedAtStart_input = payAccruedAtStart,
                         PACKAGE = "creditr")                       
  }
  return(spread)
}
