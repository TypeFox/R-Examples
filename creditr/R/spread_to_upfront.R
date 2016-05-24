#' Calculate Upfront Payments
#' 
#' \code{spread_to_upfront} takes a dataframe of variables on CDSs to return 
#' a vector of upfront values. Note that all CDS in the data frame must be denominated in 
#' the same currency.
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
#' @param notional is the amount of the underlying asset on which the
#'        payments are based. Default is 10000000, i.e. 10MM.
#' @param isPriceClean refers to the type of upfront calculated. It is
#'        boolean. When \code{TRUE}, calculate principal only. When
#'        \code{FALSE}, calculate principal + accrual.
#'        
#' @return vector of upfront values (with accrual) in the same order

spread_to_upfront <- function(x, 
                              currency.var = "currency", 
                              notional     = 10000000,
                              date.var     = "date", 
                              spread.var   = "spread",
                              coupon.var   = "coupon",
                              tenor.var    = "tenor",
                              maturity.var = "maturity",
                              recovery.var = "recovery",
                              isPriceClean = FALSE){
  
  has.maturity.var <- !(is.null(maturity.var) || is.null(x[[maturity.var]]))
  has.tenor.var <- !(is.null(tenor.var) || is.null(x[[tenor.var]]))
  stopifnot(xor(has.maturity.var, has.tenor.var))
  
  ## Coerce to desired classes
  
  x[[currency.var]] <- as.character(x[[currency.var]])
  x[[date.var]] <- as.Date(x[[date.var]])
  
  ## stop if not all the relevant variables are contained in x
  
  stopifnot(all(c(date.var, spread.var, coupon.var) %in% names(x)))
  
  ## stop if the relevant variables are not of the appropriate type 
  
  stopifnot(!has.maturity.var || inherits(x[[maturity.var]], "Date"))
  stopifnot(is.character(x[[currency.var]]))
  stopifnot(is.numeric(notional))
  stopifnot(is.numeric(x[[coupon.var]]))
  
  results <- rep(NA, nrow(x))
  
  x <- add_dates(x,
                 date.var = date.var,
                 maturity.var = maturity.var,
                 tenor.var = tenor.var,
                 currency.var = currency.var)
  
  x <- add_conventions(x)
  
  for(i in 1:nrow(x)){
    
    rates.info <- get_rates(date = as.Date(x[i, date.var]), currency = x[i, currency.var])
      
    
    
    results[i] <- .Call('calcUpfrontTest',
                        baseDate_input = separate_YMD(x$baseDate[i]),
                        types = paste(as.character(rates.info$type), collapse = ""),
                        rates = as.numeric(as.character(rates.info$rate)),
                        expiries = as.character(rates.info$expiry),
                        
                        mmDCC = x$mmDCC[i],
                        fixedSwapFreq = x$fixedFreq[i],
                        floatSwapFreq = x$floatFreq[i],
                        fixedSwapDCC = x$fixedDCC[i],
                        floatSwapDCC = x$floatDCC[i] ,
                        badDayConvZC = x$badDayConvention[i],
                        holidays = "NONE",
                        
                        todayDate_input = separate_YMD(as.Date(x[i, date.var])),
                        valueDate_input = separate_YMD(x$valueDate[i]),
                        benchmarkDate_input = separate_YMD(x$startDate[i]),
                        startDate_input = separate_YMD(x$startDate[i]),
                        endDate_input = separate_YMD(x$endDate[i]),
                        stepinDate_input = separate_YMD(x$stepinDate[i]),
                        
                        dccCDS = "ACT/360",
                        ivlCDS = "1Q",
                        stubCDS = "F",
                        badDayConvCDS = "F",
                        calendar = "NONE",
                        
                        parSpread = x[i, spread.var],
                        couponRate = x[i, coupon.var],
                        recoveryRate = x[i, recovery.var],
                        isPriceClean_input = isPriceClean,
                        payAccruedOnDefault_input = TRUE,
                        notional = notional,
                        PACKAGE = "creditr")
  } 
  return(results)
}