#### position is between 0 and 1, indicating the relative position of the ####
#### current order in the queue. ####

GetQueueLength <- function(dir = dir, date, orderTime, ticker, limitPrice, 
                           CALL, position = 1, market = 'SHSZ') {
  ## formalize argument ##
  orderTime <- as.numeric(orderTime)
  if (CALL != 'BUY' & CALL != 'SELL') {
    stop ('Wrong CALL.')
  }
  if (position < 0 | position > 1) {
    stop ('Wrong position.')
  }
  
  ## load whole day's data ##
  rawData <- LoadTickData(dir = dir, ticker = ticker, date = date, 
                          CALL = CALL, market = market)
  if (class(rawData) == 'logical') {
    return (0)
  } 
  
  ## find needed time slots ##
  data <- DataSlice(rawData, time1 = 91500, time2 = orderTime)
  if (class(data) == 'logical') {
    return (0)
  }
  
  ## exmine the bid / ask queue ##
  idx = which(data[nrow(data), 2:6] == limitPrice) 
  if (length(idx) == 0) {
    return (0)
  } else { ## idx from 1, but prices from column 2
    return (floor(data[nrow(data), idx + 6] * position / 100) * 100)
  }
 
}