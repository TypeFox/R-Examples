LoadTickData <- function(dir, ticker, date, CALL = 'BUY', market = 'SHSZ') {
  if (market == 'SHSZ') {
    return (.LoadTickDataSHSZ(dir = dir, ticker = ticker, date = date, 
                              CALL = CALL))
  } else if (market == 'HK') {
    return (.LoadTickDataHK(dir = dir, ticker = ticker, date = date, 
                            CALL = CALL))
  } else {
    return (NA)
  }
}