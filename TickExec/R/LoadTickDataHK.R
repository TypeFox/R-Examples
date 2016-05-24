#### Load Intraday Tick data, can't do time slices ####
#### CALL = 'BUY', 'SELL' or 'TRADE' ####
#### return either a data.frame, or NA ####

.LoadTickDataHK <- function(dir, ticker, date, CALL = 'BUY') {
  ## formalize arguments ##
  ticker <- as.numeric(ticker)
  ticker <- sprintf('%05d', ticker)
  date   <- as.character(date)
  CALL   <- toupper(CALL)
  
  ## formalize dir ##
  dirLast <- substr(dir, nchar(dir), nchar(dir))
  if (dirLast != '\\' & dirLast != '/') {
    dir <- paste(dir, '/', sep = '')
  }
  
  ## extract relevant columns ##
  if (CALL == 'TRADE') {
    targetFile <- paste(dir, ticker, '/', date, '/trade.csv', sep = '')
    if (file.exists(targetFile) == FALSE) { 
      return (NA) 
    } else if (file.info(targetFile)$size <= 100) {
      return (NA)
    }
    raw <- read.csv(targetFile, header = FALSE)
    raw <- raw[which(raw[, 9] != 'R'), ]
    out <- data.frame(TIME         = raw[, 3], 
                      TRADE_PRICE  = raw[, 5], 
                      TRADE_VOLUME = raw[, 7])
    out$TRADE_PRICE  <- PriceToNA(out$TRADE_PRICE)
    out$TRADE_VOLUME <- VolumeToZero(out$TRADE_VOLUME)
  } else if (CALL == 'BUY') {
    ## load bid price ##
    bidPriceFile <- paste(dir, ticker, '/', date, '/CB.csv', sep = '')
    if (file.exists(bidPriceFile) == FALSE) { 
      return (NA) 
    } else if (file.info(bidPriceFile)$size <= 100) {
      return (NA)
    }
    
    ## load bid quantity ##
    bidVolFile <- paste(dir, ticker, '/', date, '/X1.csv', sep = '')
    if (file.exists(bidVolFile) == FALSE) { 
      return (NA) 
    } else if (file.info(bidVolFile)$size <= 100) {
      return (NA)
    }
    
    bidPrice <- read.csv(bidPriceFile, header = FALSE)
    bidVol <- read.csv(bidVolFile, header = FALSE)
    
    ## extract time, price and volume## 
    time      <- bidPrice[, 1] %/% 1000
    buy1Price <- PriceToNA(bidPrice[, 3])
    buy1Vol   <- VolumeToZero(bidVol[, 3])
    
    if (length(buy1Price) != length(buy1Vol)) {
      return (NA)
    }
    
    fillPrice <- matrix(ncol = 4, nrow = length(buy1Price), NA)
    priceMatrix <- data.frame(cbind(buy1Price, fillPrice))
    colnames(priceMatrix) <- paste('BUY', 1:5, '_PRICE',    sep = '')
    
    fillVol <- matrix(ncol = 4, nrow = length(buy1Vol), 0)
    volMatrix <- data.frame(cbind(buy1Vol, fillVol))
    colnames(volMatrix) <- paste('BUY', 1:5, '_QUANTITY',    sep = '')
    
    out <- data.frame(TIME = time,
                      priceMatrix,
                      volMatrix)
  } else if (CALL == 'SELL') {
    ## load ask price ##
    askPriceFile <- paste(dir, ticker, '/', date, '/CA.csv', sep = '')
    if (file.exists(askPriceFile) == FALSE) { 
      return (NA) 
    } else if (file.info(askPriceFile)$size <= 100) {
      return (NA)
    }
    
    ## load ask quantity ##
    askVolFile <- paste(dir, ticker, '/', date, '/Y1.csv', sep = '')
    if (file.exists(askVolFile) == FALSE) { 
      return (NA) 
    } else if (file.info(askVolFile)$size <= 100) {
      return (NA)
    }
    
    askPrice <- read.csv(askPriceFile, header = FALSE)
    askVol <- read.csv(askVolFile, header = FALSE)
    
    ## extract time, price and volume## 
    time      <- askPrice[, 1] %/% 1000
    ask1Price <- PriceToNA(askPrice[, 3])
    ask1Vol   <- VolumeToZero(askVol[, 3])
    
    if (length(ask1Price) != length(ask1Vol)) {
      return (NA)
    }
    
    fillPrice <- matrix(ncol = 4, nrow = length(ask1Price), NA)
    priceMatrix <- data.frame(cbind(ask1Price, fillPrice))
    colnames(priceMatrix) <- paste('SELL', 1:5, '_PRICE',    sep = '')
    
    fillVol <- matrix(ncol = 4, nrow = length(ask1Vol), 0)
    volMatrix <- data.frame(cbind(ask1Vol, fillVol))
    colnames(volMatrix) <- paste('SELL', 1:5, '_QUANTITY',    sep = '')
    
    out <- data.frame(TIME = time,
                      priceMatrix,
                      volMatrix)
  } else {
    return (NA)
  }

  return(out)
}








