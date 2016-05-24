

PortfolioWorth <- function (dir = dir, df, date, time = 145900, market = 'SHSZ') {
  ## normalize argument ##
  if (class(df) != 'data.frame') {
    stop ('Need dataframe as input.')
  }
  if (nrow(df) == 0) {
    stop('Empty input.')
  }
  
  worth = as.numeric(sum(df$CAPLEFT) + sum(df$QUANTOUT))
  
  ## discount the holding stocks ##
  idx <- which(df$VOLUMEHOLD > 0)
  
  if (length(idx) == 0) {
    return (worth)
  }
  
  for (i in idx) {
    tempPrice = GetLastPrice(dir = dir, date = date, time = time, 
                             ticker = df[i, 'TICKER'], market = market)
    if (is.na(tempPrice) == TRUE) {
      tempPrice = df[i, 'AVGPRICEIN']
    }
    worth = worth + as.numeric(tempPrice * df[i, 'VOLUMEHOLD'])
  }
  
  return (worth)
}