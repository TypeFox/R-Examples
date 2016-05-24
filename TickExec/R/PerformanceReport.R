PerformanceReport <- function (df, cumPnL, initCap = NA) {
  ## normalize argument ##
  if (class(df) != 'data.frame') {
    stop ('Wrong data type given.')
  }
  
  if (any(is.na(cumPnL))) {
    stop ('Illegal input in PnL process.')
  }

  ## all statistic come from these two input ##
  cumPnL <- as.numeric(cumPnL)
  if (any(df$VOLUMEIN == 0)) {
    df <- df[-(which(df$VOLUMEIN == 0)), ]
  }
  if (any(df$VOLUMEHOLD > 0)) {
    df <- df[-(which(df$VOLUMEHOLD > 0)), ]
  }

  tradeDays = length(unique(df$DATEIN))
  totalDays = length(cumPnL)
  if (is.na(initCap) == TRUE) {
    dfDay1 <- df[which(df$DATEIN == df$DATEIN[1]), ]
    initCap = sum(dfDay1$QUANTIN) + sum(dfDay1$CAPLEFT)
  }
  dailyPnL <- diff(cumPnL)
  intrinsicM <- aggregate(df[, c('QUANTIN', 'QUANTOUT')], by = df['DATEIN'], 
                          FUN = sum)
  intrisicReturns <- intrinsicM$QUANTOUT / intrinsicM$QUANTIN - 1
  ## total number of days ##
  report <- data.frame(DAYS = totalDays)
  ## first and last trade ##
  report$FIRSTTRD = df[1, 'DATEIN']
  report$LASTTRD = df[nrow(df), 'DATEIN']
  ## non-trading percentage ##
  report$NONTRDPERC = 1 - tradeDays / report$DAYS
  ## trading opportunities ##
  report$DAILYTRD = nrow(df) / report$DAYS
  ## total PnL ##
  report$TOTALPNL = cumPnL[totalDays]
  perTradeReturns <- df$AVGPRICEOUT / df$AVGPRICEIN - 1
  ## per trade return ##
  report$RETPERTRD = mean(perTradeReturns)
  ## trade-wise hit rate ##
  report$TRDHITRAT = length(which(perTradeReturns > 0)) / length(perTradeReturns)
  ## daily hit rate ##
  report$DLYHITRAT = length(which(cumPnL > 0)) / totalDays
  ## simple return ##
  report$ANNRET = (cumPnL[totalDays] / initCap) / totalDays * 250
  ## sharpe ratio ##
  report$SHARPE = mean(dailyPnL) / sd(dailyPnL) * sqrt(250)
  ## draw down ##
  report$DRAWDOWN = DrawDown(cumPnL + initCap)
  ## intrinsic returns ##
  report$INRETURN = mean(intrisicReturns) * 250
  
  return(report)
}
