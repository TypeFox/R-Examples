#### If limitPrice = NA, then use the last trade price before orderFrom ####

LimitSell <- function (dir = dir, date, dfLog, limitPrice = NA, 
                       orderFrom, orderTo = 150000, orderLast = 7 * 3600,
                       costOut = 0.001, market = 'SHSZ') {
  ## formalize argumaents ##
  orderFrom = as.numeric(orderFrom)
  orderTo   = as.numeric(orderTo)
  orderLast = as.numeric(orderLast)
  ticker    = as.numeric(dfLog$TICKER)
  
  ## load whole day's data ##
  rawBidData <- LoadTickData(dir = dir, ticker = ticker, date = date, 
                             CALL = 'BUY', market = market)
  if (class(rawBidData) == 'logical') {
    return (dfLog)
  } 
  
  rawTradeData <- LoadTickData(dir = dir, ticker = ticker, date = date, 
                               CALL = 'TRADE', market = market)
  if (class(rawTradeData) == 'logical') {
    return (dfLog)
  } 
  
  ## find needed time slots ##
  bidData <- DataSlice(rawBidData, time1 = orderFrom, time2 = orderTo, 
                       last = orderLast)
  if (class(bidData) == 'logical') {
    return (dfLog)
  }
  
  tradeData <- DataSlice(rawTradeData, time1 = orderFrom, time2 = orderTo, 
                         last = orderLast)
  if (class(tradeData) == 'logical') {
    return (dfLog)
  }
  
  ## determine limit price ##
  if (is.na(limitPrice) == TRUE) {
    limitPrice = GetLastPrice(dir = dir, date = date, time = orderFrom, 
                              ticker = ticker)
  }
  if (limitPrice <= 0) {
    return (dfLog)
  }
  
  ## queue length before current order ##
  queueLength = GetQueueLength(dir = dir, date = date, orderTime = orderFrom, 
                               ticker = ticker, limitPrice = limitPrice, 
                               CALL = 'SELL', position = 1)

  ## begin execution ##
  if (max(tradeData$TRADE_PRICE) < limitPrice) {
    return (dfLog)
  }
  
  ## volume from trade ##
  tradeIdx = which(tradeData$TRADE_PRICE >= limitPrice)
  tradePotentialVol = sum(tradeData[tradeIdx, 'TRADE_VOLUME'])
  tradeExecTime = tradeData[min(tradeIdx), 'TIME']
  
  ## volume from bid ##
  bidPrices <- bidData[, 2:6]
  bidVolume <- bidData[, 7:11]
  idx2D <- which(bidPrices >= limitPrice, arr.ind = TRUE)
  if (nrow(idx2D) == 0) {
    bidPotentialVol = 0
    bidExecTime = 150000
  } else {
    uniqueBidPrice <- as.matrix(unique(bidPrices[idx2D]))
    FTempSell <- function(price) { ## to call in apply()
      idxTemp <- which(bidPrices == price, arr.ind =TRUE)
      return (max(bidVolume[idxTemp]))
    }
    bidPotentialVol = sum(apply(uniqueBidPrice, 1, FTempSell))
    
    bidExecTime = bidData[min(idx2D[, 1]), 'TIME']
  }
  
  ## calculate number of shares executed ##
  availableVol = tradePotentialVol + bidPotentialVol - queueLength
  execVol = min(dfLog$VOLUMEHOLD, availableVol)
  if (execVol <= 0) {
    return (dfLog)
  }
  
  ## average price and quantity ##
  execQuant = execVol * limitPrice / (1 + costOut)
  
  ## update log ##
  if (dfLog$SELLATTEMPT == 0) { ## first attempt to sell ##
    dfLog$VOLUMEOUT   = execVol
    dfLog$QUANTOUT    = execQuant
  } else { ## tried to sell before ##
    dfLog$VOLUMEOUT   = as.numeric(dfLog$VOLUMEOUT) + execVol
    dfLog$QUANTOUT    = as.numeric(dfLog$QUANTOUT) + execQuant
  }
  dfLog$AVGPRICEOUT = dfLog$QUANTOUT / dfLog$VOLUMEOUT
  dfLog$DATEOUT = date
  dfLog$TIMEOUT = min(bidExecTime, tradeExecTime)
  dfLog$DEPTHOUT = 0
  dfLog$SELLATTEMPT = dfLog$SELLATTEMPT + 1
  dfLog$VOLUMEHOLD  = dfLog$VOLUMEHOLD - execVol
  
  return(dfLog)
}




