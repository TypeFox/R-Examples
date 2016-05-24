

GetLastPrice <- function(dir = dir, date, time, ticker, market = 'SHSZ') {
  ## formalize argumaents ##
  time <- as.numeric(time)
  
  ## load whole day's data ##
  rawData <- LoadTickData(dir = dir, ticker = ticker, date = date, 
                          CALL = 'TRADE', market = market)
  if (class(rawData) == 'logical') {
    return (NA)
  } 
  
  ## find needed time slots ##
  data <- DataSlice(rawData, time1 = 90000, time2 = time)
  if (class(data) == 'logical') {
    return (NA)
  }
  
  ## retrive wanted price ##
  return (data$TRADE_PRICE[nrow(data)])
}