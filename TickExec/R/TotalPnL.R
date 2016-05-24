TotalPnL <- function(dir = dir, df, date, time = 160000, market = 'SHSZ') {
  ## normalize argument ##
  if (class(df) != 'data.frame') {
    stop ('Wrong data type given.')
  }
  
  if (nrow(df) == 0) {
    return (0)
  }
  
  worth <- PortfolioWorth(dir = dir, df = df, date = date, time = time,
                          market = market)
  
  worth = worth - sum(df$CAPLEFT, na.rm = TRUE) - sum(df$QUANTIN, na.rm = TRUE)
  return(as.numeric(worth))
}