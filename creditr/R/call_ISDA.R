#' call ISDA c function
#' 
#' \code{call_ISDA} call ISDA function
#' 
#' @param x dataframe which contains relevant dates and convention info
#' @param name character function name within which call_ISDA is called
#' @param rates.info dataframe which contains relevant rates data
#'   
#' @return a numeric value which is the difference between the new upfront and 
#'   the old one

call_ISDA <- function(x, name, rates.info){

  stopifnot(is.data.frame(x))
  stopifnot(c("date", "currency", "spread", "coupon",
              "recovery", "notional", "stepinDate", "valueDate",
              "startDate", "firstcouponDate", "pencouponDate", "endDate",
              "backstopDate", "baseDate", "badDayConvention", "mmDCC",
              "mmCalendars", "fixedDCC", "floatDCC", "fixedFreq",
              "floatFreq", "swapCalendars") %in% names(x))
  
  stopifnot(is.data.frame(rates.info))
  stopifnot(c("expiry", "rate", "type") %in% names(rates.info))
  
  ## define CS10, IR.DV01, rec.risk.01 and spread.DV01 and set their default
  ## to zero. If the name matches any of the four functions, then the corres-
  ## ponding value for that function will be changed to one. As can be seen
  ## in the code below, this does not mean that the value of these inputs to C
  ## code are changed to one, but rather that either 1% is added or 10% is
  ## multiplied.
  
  CS10        <- 0
  IR.DV01     <- 0
  rec.risk.01 <- 0
  spread.DV01 <- 0
  
  if(name == "CS10"){                ## if CS10 is to be changed,
    CS10 <- 1                        ## then 
  } else if(name == "IR.DV01"){
    IR.DV01 <- 1
  } else if(name == "rec.risk.01"){
    rec.risk.01 <- 1
  } else if(name == "spread.DV01"){
    spread.DV01 <- 1
  } else{
    warning("No change provided")
  }
  
  upfront.new <- .Call('calcUpfrontTest',
                       baseDate_input = separate_YMD(x$baseDate),
                       types = paste(as.character(rates.info$type), collapse = ""),
                       
                       ## If IRDV01 is changed, then rates is changed by adding 1/1e4
                       ## otherwise it remains unchanged
                       
                       rates = as.numeric(as.character(rates.info$rate)) + IR.DV01 * 1/1e4,
                       expiries = as.character(rates.info$expiry),
                       
                       mmDCC = as.character(x$mmDCC),
                       fixedSwapFreq = as.character(x$fixedFreq),
                       floatSwapFreq = as.character(x$floatFreq),
                       fixedSwapDCC = as.character(x$fixedDCC),
                       floatSwapDCC = as.character(x$floatDCC),
                       badDayConvZC = as.character(x$badDayConvention),
                       holidays = "None",
                       
                       todayDate_input = separate_YMD(x$date),
                       valueDate_input = separate_YMD(x$valueDate),
                       benchmarkDate_input = separate_YMD(x$startDate),
                       startDate_input = separate_YMD(x$startDate),
                       endDate_input = separate_YMD(x$endDate),
                       stepinDate_input = separate_YMD(x$stepinDate),
                       
                       dccCDS = "ACT/360",
                       ivlCDS = "1Q",
                       stubCDS = "F",
                       badDayConvCDS = "F",
                       calendar = "None",
                       
                       ## If spread.DV01 remains unchanged and 
                       ##CS10 is changed, then the spread is changed 
                       ## by multipling 1.1; otherwise, spread remains 
                       ## unchanged.
                       
                       ## If CS10 remains unchanged and spread.DV01 is
                       ## changed, then the spread is changed by adding 1;
                       ## otherwise it remains unchanged.
                       
                       parSpread = x$spread * (1 + CS10 * 0.1) + spread.DV01,
                       couponRate = x$coupon,
                       
                       ## if rec.risk is changed, then recovery rate is
                       ## is changed by adding 1%. otherwise it remains
                       ## unchanged
                       
                       recoveryRate = x$recovery + 0.01 * rec.risk.01,
                       isPriceClean_input = FALSE,
                       payAccruedOnDefault_input = TRUE,
                       notional = x$notional,
                       PACKAGE = "creditr")
  
  upfront.orig <- .Call('calcUpfrontTest',
                        baseDate_input = separate_YMD(x$baseDate),
                        types = paste(as.character(rates.info$type), collapse = ""),
                        rates = as.numeric(as.character(rates.info$rate)),
                        expiries = as.character(rates.info$expiry),
                        
                        mmDCC = as.character(x$mmDCC),
                        fixedSwapFreq = as.character(x$fixedFreq),
                        floatSwapFreq = as.character(x$floatFreq),
                        fixedSwapDCC = as.character(x$fixedDCC),
                        floatSwapDCC = as.character(x$floatDCC),
                        badDayConvZC = as.character(x$badDayConvention),
                        holidays = "None",
                        
                        todayDate_input = separate_YMD(x$date),
                        valueDate_input = separate_YMD(x$valueDate),
                        benchmarkDate_input = separate_YMD(x$startDate),
                        startDate_input = separate_YMD(x$startDate),
                        endDate_input = separate_YMD(x$endDate),
                        stepinDate_input = separate_YMD(x$stepinDate),
                        
                        dccCDS = "ACT/360",
                        ivlCDS = "1Q",
                        stubCDS = "F",
                        badDayConvCDS = "F",
                        calendar = "None",
                        
                        parSpread = x$spread,
                        couponRate = x$coupon,
                        recoveryRate = x$recovery,
                        isPriceClean_input = FALSE,
                        payAccruedOnDefault_input = TRUE,
                        notional = x$notional,
                        PACKAGE = "creditr")
  
  return(upfront.new - upfront.orig)
}