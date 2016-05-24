

MarketBuy <- function (dir = dir, date, ticker, capital, orderTime, 
                       costIn = 0.001, market = 'SHSZ') {
  ## formalize argumaents ##
  orderTime <- as.numeric(orderTime)
  
  ## load whole day's data ##
  rawData <- LoadTickData(dir = dir, ticker = ticker, date = date, 
                          CALL = 'SELL', market = market)
  if (class(rawData) == 'logical') {
    return (InitLogEntry(dateIn = date, ticker = ticker, capital = capital))
  } 
  
  ## find needed time slots ##
  data <- DataSlice(rawData, time1 = orderTime, last = 150000 - orderTime)
  if (class(data) == 'logical') {
    return (InitLogEntry(dateIn = date, ticker = ticker, capital = capital))
  }
  
  ## locate the needed time slice ##
  entry <- data[1, ]
  
  ## calculate average price and depth ##
  prices <- as.numeric(entry[1, 2:6]) * (1 + costIn)
  volumes <- as.numeric(entry[1, 7:11])
  if (sum(volumes) == 0){
    return (InitLogEntry(dateIn = date, ticker = ticker, capital = capital))
  }
  
  potentialQuant    <- VolumeToZero(prices * volumes)
  cumPotentialQuant <- cumsum(potentialQuant)
  
  depth = min(which(cumPotentialQuant > capital), 6)
  if (depth == 1) {
    execVol   = floor(capital / prices[1] / 100) * 100
    execQuant = prices[1] * execVol
  } else if (depth <= 5) {
    remainderQuant = capital - cumPotentialQuant[depth - 1]
    remainderVol   = floor(remainderQuant / prices[depth] / 100) * 100
    execVol        = sum(volumes[1:(depth - 1)]) + remainderVol
    execQuant      = cumPotentialQuant[depth - 1] + remainderVol * prices[depth]
  } else if (depth == 6) {
    depth     = 5
    execQuant = sum(potentialQuant)
    execVol   = sum(volumes)
  }
  
  if (execVol == 0){
    return (InitLogEntry(dateIn = date, ticker = ticker, capital = capital))
  }
  avgPrice  = execQuant / execVol
  
  return (InitLogEntry(dateIn    = date, 
                         timeIn    = entry$TIME, 
                         ticker    = ticker, 
                         capital   = capital,
                         execVol   = execVol,
                         execQuant = execQuant,
                         avgPrice  = avgPrice,
                         depthIn   = depth))
}









