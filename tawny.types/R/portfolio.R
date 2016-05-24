##---------------------------- TAWNY_PORTFOLIO -----------------------------##
# Example
# p <- TawnyPortfolio(c('FCX','AAPL','JPM','AMZN','VMW','TLT','GLD','FXI','ILF','XOM'))
TawnyPortfolio(symbols, window, obs) %::% character : numeric : numeric : list
TawnyPortfolio(symbols, window=90, obs=150) %as%
{
  returns <- AssetReturns(symbols, obs=obs)
  TawnyPortfolio(returns, window)
}


TawnyPortfolio(returns, window, extra) %::% AssetReturns : numeric : . : list
TawnyPortfolio(returns, window=90, extra=NULL) %as%
{
  periods <- anylength(returns) - window + 1
  c(list(symbols=anynames(returns), window=window, obs=anylength(returns),
       periods=periods, returns=returns), extra)
}


TawnyPortfolio(returns, window, extra) %::% zoo : numeric : . : list
TawnyPortfolio(returns, window, extra=NULL) %as%
{
  TawnyPortfolio(AssetReturns(returns), window, extra)
}



rollapply.TawnyPortfolio <- function(x, fun, ...)
{
  steps <- array(seq(1,x$periods))
  out <- apply(steps, 1, function(idx) fun(window_at(x,idx), ...))
  out <- t(out)
  rownames(out) <- format(index(x$returns[(x$obs - x$periods + 1):x$obs,]))
  out
}

window_at(x, idx) %::% TawnyPortfolio : a : TawnyPortfolio
window_at(x, idx) %as%
{
  returns <- x$returns[idx:(x$window + idx - 1),]
  x$returns <- returns
  x
}

start.TawnyPortfolio <- function(x, ...)
{
  start(x$returns[x$obs - x$periods,])
}

end.TawnyPortfolio <- function(x, ...)
{
  end(x$returns)
}


##------------------------- BENCHMARK PORTFOLIO ----------------------------##
# Convenience function for creating a benchmark portfolio
# m <- BenchmarkPortfolio('^GSPC', 150, 200)
BenchmarkPortfolio(market,window,obs, end=Sys.Date(),...) %as%
{
  if (is.character(market))
  {
    start <- end - (10 + obs * 365/250)
    mkt <- getSymbols(market, src='yahoo',from=start,to=end, auto.assign=FALSE)
  }
  else mkt <- market
  # xts has moved to use POSIX
  #end <- as.POSIXct(end)

  mkt.ret <- Delt(Cl(mkt))
  mkt.ret <- mkt.ret[index(mkt.ret) <= end]
  mkt.ret <- tail(mkt.ret, obs)
  colnames(mkt.ret) <- 'benchmark'

  w.count <- obs - window + 1
  weights <-
    xts(matrix(1,ncol=1,nrow=w.count), order.by=index(tail(mkt.ret, w.count)))

  TawnyPortfolio(mkt.ret, window, list(rf.rate=0.01))
}


# Calculate portfolio returns based
# returns <- PortfolioReturns(p, weights)
# chart.PerformanceSummary(returns)
PortfolioReturns(p, weights) %::% TawnyPortfolio : numeric : a
PortfolioReturns(p, weights) %as%
{
  PortfolioReturns(p$returns, weights)
}

PortfolioReturns(h, weights) %::% AssetReturns : numeric : a
PortfolioReturns(h, weights) %as%
{
  # Shift dates so weights are used on following date's data for out-of-sample
  # performance
  w.index <- c(index(weights[2:anylength(weights)]), end(weights) + 1)
  index(weights) <- w.index
 
  h.trim <- h[index(h) %in% index(weights)]
  ts.rets <- apply(zoo(h.trim) * zoo(weights), 1, sum)

  # This is in here to fix some strange behavior related to rownames vs index
  # in zoo objects and how they are used after an apply function
  #names(ts.rets) <- index(h.trim)
  #ts.rets <- zoo(ts.rets, order.by=as.Date(names(ts.rets)))
  ts.rets <- zoo(ts.rets, order.by=index(h.trim))
  
  # This causes problems
  if (any(is.na(ts.rets)))
  {
    flog.warn("Filling NA returns with 0")
    ts.rets[is.na(ts.rets)] <- 0
  }
  
  flog.debug("Returns count: %s", anylength(ts.rets))

  return(ts.rets)
}


# This produces a portfolio in matrix format (t x m) as a zoo class. 
# Params
#  symbols: A vector of symbols to retrieve. This uses quantmod to retrieve
#    the data.
#  obs: The number of observations that you want. Use this if you want the 
#    number of points to be explicit. Either obs or start is required.
#  start: The start date, if you know that explicitly. Using this will ensure
#    that the data points are bound to the given range but the precise number
#    of points will be determined by the number of trading days.
#  end: The most recent date of observation. Defaults to current day.
#  fun: A function to use on each symbol time series. Defaults to Cl to operate
#    on close data. For expected behavior, your function should only return
#    one time series.
# TODO: 
#  Fix names
#  Add method to add other portfolio elements (such as synthetic securities)
# Example:
#  h <- AssetReturns(c('GOOG','AAPL','BAC','C','F','T'), 150)
AssetReturns(returns) %::% zoo : zoo
AssetReturns(returns) %as% returns

AssetReturns(symbols, obs=NULL, start=NULL, end=Sys.Date(),
  fun=function(x) Delt(Cl(x)), reload=FALSE, na.value=NA, ...) %as%
{
  if (is.null(start) & is.null(obs)) { stop("Either obs or start must be set") }
  end <- as.Date(end)

  # Estimate calendar days from windowed business days. The 10 is there to
  # ensure enough points, which get trimmed later
  if (is.null(start)) { start <- end - (10 + obs * 365/250) }

  #ensure(symbols, src='yahoo', reload=reload, from=start, to=end, ...)

  # Merge into a single zoo object
  p <- xts(order.by=end)
  for (s in symbols)
  {
    asset <- getSymbols(s, from=start, to=end, auto.assign=FALSE)
    raw <- fun(asset)
    flog.info("Binding %s for [%s,%s]",s, format(start(raw)),format(end(raw)))
      
    a <- xts(raw, order.by=index(asset))
    p <- cbind(p, a[2:anylength(a)])
  }
  colnames(p) <- symbols
  # First remove dates that have primarily NAs (probably bad data)
  o.dates <- rownames(p)
  p <- p[apply(p, 1, function(x) sum(x, na.rm=TRUE) != 0), ]
  flog.info("Removed suspected bad dates %s",setdiff(o.dates,rownames(p)))

  if (! is.na(na.value))
  {
    #for (s in symbols) p[,s][is.na(p[,s])] <- na.value
    p[is.na(p)] <- 0
    flog.info("Replaced NAs with %s",na.value)
  }
  else
  {
    # NOTE: This has consistency issues when comparing with a market index
    o.dates <- rownames(p)
    p <- p[apply(p, 1, function(x) sum(is.na(x)) < 0.1 * length(x) ), ]
    flog.info("Removed dates with too many NAs %s",setdiff(o.dates,rownames(p)))

    # Now remove columns with NAs
    nas <- apply(p, 2, function(x) !any(is.na(x)) )
    p <- p[,which(nas == TRUE)]
    flog.info("Removed symbols with NAs: %s",setdiff(symbols,anynames(p)))
  }

  if (is.null(obs)) { return(p[paste(start,end, sep='::')]) }

  p <- p[index(p) <= end]
  idx.inf <- anylength(p) - min(anylength(p), obs) + 1
  idx.sup <- anylength(p)
  
  flog.info("Loaded portfolio with %s assets",ncol(p))
  out <- p[idx.inf:idx.sup, ]
  class(out) <- c('returns', class(out))

  if (is.null(rownames(out))) rownames(out) <- format(index(out), "%Y-%m-%d")
  out
}



# Generate the composition for an equity index
# Example
# Get SP500 components
#   sp500.idx <- EquityIndex()
# Get DOW components
#   dow.idx <- EquityIndex('^DJI')
# Get FTSE components
#   ftse.idx <- EquityIndex('^FTSE')
# Get HSI components
#   hsi.idx <- EquityIndex('^HSI')
# h <- AssetReturns(EquityIndex('^DJI'), obs=100)
EquityIndex(ticker='^GSPC', hint=NA, src='yahoo') %as%
{
  if (is.na(hint))
  {
    hints <- c(500, 30, 102, 42)
    names(hints) <- c('^GSPC', '^DJI', '^FTSE', '^HSI')
    hint <- hints[ticker]
  }

  # http://download.finance.yahoo.com/d/quotes.csv?s=@%5EGSPC&f=sl1d1t1c1ohgv&e=.csv&h=0
  # TODO Fix the URL for the composition
  base <- 'http://download.finance.yahoo.com/d/quotes.csv?s=@'
  formats <- '&f=sl1d1t1c1ohgv&e=.csv&h='

  comp <- NULL
  pages = max(1, hint %/% 50)
  for (page in 1:pages)
  {
    start <- (page-1) * 50 + 1
    url <- paste(base, ticker, formats, start, sep='')
    flog.info("Loading page %s for %s",page,ticker)
    data <- read.csv(url, header=FALSE)

    # This is here due to a bug in Yahoo's download where the first record gets
    # duplicated in each subsequent page
    idx = 2; if (page == 1) { idx = 1 }
    comp <- rbind(comp, data[idx:anylength(data),])

  }
  as.character(comp[,1])
}


