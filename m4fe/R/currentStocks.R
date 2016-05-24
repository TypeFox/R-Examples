#' @title Current Stock Info (Yahoo Finance API)
#' @description Parses the Yahoo Finance API to retrieve current stock information
#' @param stockIndex The name of the stock index to retrieve information for (default GOOGL)
#' @details Uses R's read.csv method and generates a data table containing stock information such as Bid-Ask, EPS.
#' @examples currentStocks()
#' @export
currentStocks <- function(stockIndex="GOOGL") {
  # OPTION PARAMETERS FOR YAHOO API
  pricingOptions = c('a','b','b2','b3','p','o')
  dividendOptions = c('y','d','r1','q')
  dateOptions = c('c1','c','c6','k2','p2','d1','d2','t1')
  averagesOptions = c('c8','c3','m5','m6','m7','m8','m3','m4','g','h','k1','l','l1', 't8')
  miscOptions = c('w1','w4','p1','m','m2','g1','g3','g4','g5','g6','t7','t6','i5','l2','l3','v1','v7','s6')
  weekOptions = c('k','j','j5','k4','j6','k5','w')
  symbolOptions = c('v','j1','j3','f6','n','n4','s','s1','x','s2')
  ratioOptions = c('e','e7','e8','e9','b4','j4','p5','p6','r','r2','r5','r6','r7','s7')
  
  # DESCRIPTIONS FOR OPTIONS PARAMETERS
  pricingOptionsNames = c('ask', 'bid', 'ask (realtime)', 'bid (realtime)', 'previous close', 'open')
  dividendOptionsNames = c('dividend yield', 'dividend per share', 'dividend pay rate', 'ex-dividend rate')
  dateOptionsNames = c('change', 'change and % change', 'change (realtime)', 'change percent (realtime)', 'change in percent', 'last trade date', 'trade date', 'last trade time')
  averagesOptionsNames = c('after hours change (realtime)', 'commission', 'change from 200 day moving avg', '% change from 200 day moving avg', 'change from 50 day moving avg', '% change from 50 day moving avg', '50 day moving day average', '200 day moving average', 'Todays low', 'Todays high', 'Last trade (realtime) with time', 'Last trade (with time)', 'Last trade (price only)', '1-year target price')
  miscOptionsNames = c('Days value change', 'Days value change (realtime)', 'price paid', 'days range', 'days range (realtime)', 'holdings gain percent', 'annualized gain', 'holdings gain', 'holdings gain percent (realtime)', 'holding gain (realtime)', 'ticker trend', 'trade links', 'order book (realtime)', 'high limit', 'low limit', 'holdings value', 'holdings value (realtime)', 'revenue')
  weekOptionsNames = c('52week high', '52week low', 'change from 52week low', 'change from 52week high', '% change from 52week low', '% change from 52week high', '52week range')
  symbolOptionsNames = c('more info', 'market capitalization', 'market cap (realtime)', 'float shares', 'name', 'notes', 'symbol', 'shares owned', 'stock exchange', 'shares outstanding')
  ratioOptionsNames = c('eps', 'eps estimate current year', 'eps estimate next year', 'eps estimate next quarter', 'book-value', 'EBITDA', 'price/sales', 'price/book', 'p/e ratio', 'p/e ratio (realtime', 'peg ratio', 'price/eps estimate current year', 'price/eps estimate next year', 'short ratio')
  
  # COMBINE OPTION AND DESCRIPTION ARRAYS
  options=c(pricingOptions,dividendOptions,dateOptions,averagesOptions,miscOptions,weekOptions,symbolOptions,ratioOptions)
  names=c(pricingOptionsNames,dividendOptionsNames,dateOptionsNames,averagesOptionsNames,miscOptionsNames,weekOptionsNames,symbolOptionsNames,ratioOptionsNames)
  
  rm(pricingOptions,dividendOptions,dateOptions,averagesOptions,miscOptions,weekOptions,symbolOptions,ratioOptions)
  rm(pricingOptionsNames,dividendOptionsNames,dateOptionsNames,averagesOptionsNames,miscOptionsNames,weekOptionsNames,symbolOptionsNames,ratioOptionsNames)
  
  # CREATE THE URL
  base = "http://finance.yahoo.com/d/quotes.csv"
  symbols = paste("?s=", stockIndex, sep="")
  info = "&f="
  
  for (i in 1:length(options)) info = paste(info,options[i],sep="")
  url = paste(base,symbols,info,sep="")
  
  # PARSE THE CSV FILE
  data = read.csv(url, header=FALSE)
  
  stockInfo <- character(0)
  for(i in 1:length(data)) { 
    stockInfo = c(stockInfo, data[1,i])
  }
  rm(data,base,symbols,info,i)
  
  index=data.frame(names,stockInfo[1:length(names)])
  names(index)=c('category','value')
  index
}