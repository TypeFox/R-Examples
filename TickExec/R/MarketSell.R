

MarketSell <- function(dir = dir, date, orderTime, dfLog, costOut = 0.001,
                       market = 'SHSZ') {
  ## formalize argumaents ##
  orderTime <- as.numeric(orderTime)
  ticker = as.numeric(dfLog$TICKER)
  if (dfLog$VOLUMEHOLD <= 0) {
    return (dfLog)
  } 
  
  ## load whole day's data ##
  rawData <- LoadTickData(dir = dir, ticker = ticker, date = date, 
                          CALL = 'BUY', market = market)
  if (class(rawData) == 'logical') {
    return (dfLog)
  } 
  
  ## find needed time slots ##
  data <- DataSlice(rawData, time1 = orderTime, last = 150000 - orderTime)
  if (class(data) == 'logical') {
    return (dfLog)
  }
  
  ## locate the needed time slice ##
  entry <- data[1, ]
  
  ## calculate average price and depth ##
  prices <- as.numeric(entry[1, 2:6]) / (1 + costOut)
  volumes <- as.numeric(entry[1, 7:11])
  if (sum(volumes) == 0){
    return (dfLog)
  }
  
  ## temporary summary statistics ##
  potentialQuant <- VolumeToZero(prices * volumes)
  cumPotentialVol <- cumsum(volumes)
  volToSell = as.numeric(dfLog$VOLUMEHOLD)
  depth = min(which(cumPotentialVol > volToSell), 6)
  
  if (dfLog$SELLATTEMPT == 0) { ## first attempt to sell ##
    if (depth == 1) {
      dfLog$AVGPRICEOUT = prices[1]
      dfLog$VOLUMEOUT   = volToSell
      dfLog$DEPTHOUT    = 1
      dfLog$VOLUMEHOLD  = 0
    } else if (depth <= 5) {
      remainderVol = volToSell - sum(volumes[1:(depth - 1)])
      quantOut = sum(potentialQuant[1:(depth - 1)]) + 
                 prices[depth] * remainderVol
      dfLog$AVGPRICEOUT = quantOut / volToSell
      dfLog$VOLUMEOUT   = volToSell
      dfLog$DEPTHOUT    = depth
      dfLog$VOLUMEHOLD  = 0
    } else if (depth == 6) {
      dfLog$AVGPRICEOUT = sum(potentialQuant) / sum(volumes)
      dfLog$VOLUMEOUT   = sum(volumes)
      dfLog$DEPTHOUT    = 5
      dfLog$VOLUMEHOLD  = volToSell - dfLog$VOLUMEOUT
    }
  } else { ## tried to sell before ##
    dfTemp <- dfLog[, c('VOLUMEOUT', 'QUANTOUT')]
    if (depth == 1) {
      dfLog$VOLUMEOUT   = as.numeric(dfTemp$VOLUMEOUT) + volToSell
      dfLog$AVGPRICEOUT = (dfTemp$QUANTOUT + volToSell * prices[1]) / 
                          dfLog$VOLUMEOUT 
      dfLog$DEPTHOUT    = 1
      dfLog$VOLUMEHOLD  = 0
    } else if (depth <= 5) {
      remainderVol = volToSell - sum(volumes[1:(depth - 1)])
      quantOut = sum(potentialQuant[1:(depth - 1)]) + 
        prices[depth] * remainderVol
      dfLog$VOLUMEOUT   = as.numeric(dfTemp$VOLUMEOUT) + volToSell
      dfLog$AVGPRICEOUT = (quantOut + dfTemp$QUANTOUT) / dfLog$VOLUMEOUT
      dfLog$DEPTHOUT    = depth
      dfLog$VOLUMEHOLD  = 0
    } else if (depth == 6) {
      dfLog$VOLUMEOUT   = as.numeric(dfTemp$VOLUMEOUT) + sum(volumes)
      dfLog$AVGPRICEOUT = (dfTemp$QUANTOUT + sum(potentialQuant)) / 
                          dfLog$VOLUMEOUT
      dfLog$DEPTHOUT    = 5
      dfLog$VOLUMEHOLD  = volToSell - sum(volumes)
    }
  }
  dfLog$DATEOUT = date
  dfLog$TIMEOUT = entry$TIME
  dfLog$QUANTOUT = dfLog$AVGPRICEOUT * dfLog$VOLUMEOUT
  dfLog$SELLATTEMPT = dfLog$SELLATTEMPT + 1
  return(dfLog)
}