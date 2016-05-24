#### Load Intraday Tick data, can't do time slices ####
#### CALL = 'BUY', 'SELL' or 'TRADE' ####
#### return either a data.frame, or NA ####

.LoadTickDataSHSZ <- function(dir, ticker, date, CALL = 'BUY') {
  ## formalize arguments ##
  ticker <- as.numeric(ticker)
  ticker <- sprintf('%06d', ticker)
  date   <- as.character(date)
  CALL   <- toupper(CALL)
  
  ## formalize dir ##
  dirLast <- substr(dir, nchar(dir), nchar(dir))
  if (dirLast != '\\' & dirLast != '/') {
    dir <- paste(dir, '/', sep = '')
  }
  
  ## load raw data ##
  targetFile <- paste(dir, ticker, '/', date, '.csv', sep = '')
  if (file.exists(targetFile) == FALSE) { 
    return (NA) 
  } else if (file.info(targetFile)$size <= 100) {
    return (NA)
  }
  raw <- read.csv(targetFile, header = TRUE)
  colnames(raw) <- toupper(colnames(raw))
  
  ## extract relevant columns ##
  if (CALL == 'TRADE') {
    out <- data.frame(TIME         = raw$TIME, 
                      TRADE_PRICE  = raw$TRADE_PRICE, 
                      TRADE_VOLUME = raw$TRADE_VOLUME * 100)
    out$TRADE_PRICE  <- PriceToNA(out$TRADE_PRICE)
    out$TRADE_VOLUME <- VolumeToZero(out$TRADE_VOLUME)
  } else {
    relatedPrice    <- paste(CALL, 1:5, '_PRICE',    sep = '')
    relatedQuantity <- paste(CALL, 1:5, '_QUANTITY', sep = '')
    out <- data.frame(TIME = raw$TIME, 
                      raw[ ,relatedPrice], 
                      raw[, relatedQuantity] * 100)
    out[, relatedPrice]    <- PriceToNA(out[, relatedPrice])
    out[, relatedQuantity] <- VolumeToZero(out[, relatedQuantity])
  }
  return(out)
}