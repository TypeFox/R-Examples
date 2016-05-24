#### If limitPrice = NA, then use the last trade price before orderFrom ####

LimitBuy <- function (dir = dir, date, ticker, capital, limitPrice = NA, 
                      orderFrom, orderTo = 150000, orderLast = 7 * 3600,
                      costIn = 0.001, market = 'SHSZ') {
  ## formalize argumaents ##
  orderFrom = as.numeric(orderFrom)
  orderTo   = as.numeric(orderTo)
  orderLast = as.numeric(orderLast)
  
  ## load whole day's data ##
  rawAskData <- LoadTickData(dir = dir, ticker = ticker, date = date, 
                             CALL = 'SELL', market = market)
  if (class(rawAskData) == 'logical') {
    return (InitLogEntry(dateIn = date, ticker = ticker, capital = capital))
  } 
  
  rawTradeData <- LoadTickData(dir = dir, ticker = ticker, date = date, 
                               CALL = 'TRADE', market = market)
  if (class(rawTradeData) == 'logical') {
    return (InitLogEntry(dateIn = date, ticker = ticker, capital = capital))
  } 
  
  ## find needed time slots ##
  askData <- DataSlice(rawAskData, time1 = orderFrom, time2 = orderTo, 
                       last = orderLast)
  if (class(askData) == 'logical') {
    return (InitLogEntry(dateIn = date, ticker = ticker, capital = capital))
  }
  
  tradeData <- DataSlice(rawTradeData, time1 = orderFrom, time2 = orderTo, 
                         last = orderLast)
  if (class(tradeData) == 'logical') {
    return (InitLogEntry(dateIn = date, ticker = ticker, capital = capital))
  }
  
  ## determine limit price ##
  if (is.na(limitPrice) == TRUE) {
    limitPrice = GetLastPrice(dir = dir, date = date, time = orderFrom, 
                              ticker = ticker)
  }
  if (limitPrice <= 0) {
    return (InitLogEntry(dateIn = date, ticker = ticker, capital = capital))
  }
  
  ## queue length before current order ##
  queueLength = GetQueueLength(dir = dir, date = date, orderTime = orderFrom, 
                               ticker = ticker, limitPrice = limitPrice, 
                               CALL = 'BUY', position = 1)
  
  ## begin execution ##
  if (min(tradeData$TRADE_PRICE) > limitPrice) {
    return (InitLogEntry(dateIn = date, ticker = ticker, capital = capital))
  }
  
  ## volume from trade ##
  tradeIdx = which(tradeData$TRADE_PRICE <= limitPrice)
  tradePotentialVol = sum(tradeData[tradeIdx, 'TRADE_VOLUME'])
  tradeExecTime = tradeData[min(tradeIdx), 'TIME']
    
  ## volume from ask ##
  askPrices <- askData[, 2:6]
  askVolume <- askData[, 7:11]
  idx2D <- which(askPrices <= limitPrice, arr.ind = TRUE)
  if (nrow(idx2D) == 0) {
    askPotentialVol = 0
    askExecTime = 150000
  } else {
    uniqueAskPrice <- as.matrix(unique(askPrices[idx2D]))
    FTempBuy <- function(price) { ## to call in apply()
      idxTemp <- which(askPrices == price, arr.ind =TRUE)
      return (max(askVolume[idxTemp]))
    }
    askPotentialVol = sum(apply(uniqueAskPrice, 1, FTempBuy))
    
    askExecTime = askData[min(idx2D[, 1]), 'TIME']
  }
  
  ## calculate number of shares executed ##
  idealVol = floor(capital / limitPrice / (1 + costIn) / 100) * 100
  availableVol = tradePotentialVol + askPotentialVol - queueLength
  execVol = min(idealVol, availableVol)
  if (execVol <= 0) {
    return (InitLogEntry(dateIn = date, ticker = ticker, capital = capital))
  }
  
  ## average price and quantity ##
  avgPrice = limitPrice * (1 + costIn)
  execQuant = execVol * avgPrice
  execTime = min(askExecTime, tradeExecTime)
  
  return (InitLogEntry(dateIn    = date, 
                       timeIn    = execTime, 
                       ticker    = ticker, 
                       capital   = capital,
                       execVol   = execVol,
                       execQuant = execQuant,
                       avgPrice  = avgPrice,
                       depthIn   = 0))
}



