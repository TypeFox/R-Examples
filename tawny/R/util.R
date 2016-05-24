# Look at rmetrics.org
#library(PerformanceAnalytics)
#library(futile)
#library(zoo)
#library(quantmod)

# Assign a class to the given object
# If h is not square, then assume it is returns
# If h has values > 1, assume covariance matrix of returns
# Else assume correlation matrix of returns
# OBSOLETE in 2.0
#classify <- function(x)
#{
#  if (is.null(dim(x))) stop("h must have a dim attribute")
#
#  tawny.types <- c('returns','covariance','correlation')
#  if (any(tawny.types %in% class(x))) invisible(x)
#
#  if (ncol(x) != nrow(x)) class(x) <- c(class(x), 'returns')
#  else if (max(x) > 1) class(x) <- c(class(x), 'covariance')
#  else class(x) <- c(class(x), 'correlation')
#  invisible(x)
#}



##---------------------------- PUBLIC FUNCTIONS -----------------------------##
# Optimize a portfolio using the specified portfolio generator and correlation
# matrix generator. This is a generic optimizer that allows any custom generator
# of correlation matrices to be used.


# Ensures that a given series exists and downloads via quantmod if it doesn't
# serie - Either a string or a list
ensure <- function(serie, src='FRED', reload=FALSE, ...)
{
  for (series in serie)
  {
    if (exists('from')) from <- get('from')
    if (exists('to')) to <- get('to')

    # Force a reload under certain situations
    cleaned <- sub('^','', series, fixed=TRUE)
    if (! exists(cleaned)) reload <- TRUE
    else if (! 'zoo' %in% class(get(cleaned)) ) reload <- TRUE
    else if (! 'Date' %in% class(start(get(cleaned))) ) reload <- TRUE
    else if (! 'Date' %in% class(end(get(cleaned))) ) reload <- TRUE
    else if (exists('from')) { if (start(get(cleaned)) > from) reload <- TRUE }
    else if (exists('to')) { if (end(get(cleaned)) < to) reload <- TRUE }

    if (! reload) next

    flog.debug("(Re)loading symbol %s from %s",series,src)
    getSymbols(series, src=src, ...)
  }
}

# Example
# Get SP500 components
#   sp500.idx <- getIndexComposition()
# Get DOW components
#   dow.idx <- getIndexComposition('^DJI')
# Get FTSE components
#   ftse.idx <- getIndexComposition('^FTSE')
# Get HSI components
#   hsi.idx <- getIndexComposition('^HSI')
# h <- getPortfolioReturns(getIndexComposition('^DJI'), obs=100)
getIndexComposition <- function(ticker='^GSPC', hint=NA, src='yahoo')
{
  if (is.na(hint))
  {
    hints <- c(500, 30, 102, 42)
    names(hints) <- c('^GSPC', '^DJI', '^FTSE', '^HSI')
    hint <- hints[ticker]
  }

  # http://download.finance.yahoo.com/d/quotes.csv?s=@%5EGSPC&f=sl1d1t1c1ohgv&e=.csv&h=0
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
getPortfolioReturns <- function(symbols, obs=NULL, start=NULL, end=Sys.Date(),
  fun=function(x) Delt(Cl(x)), reload=FALSE, na.value=NA, ...)
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
  #if (logLevel() > 0) cat("Returning rows [",idx.inf,",",idx.sup,"]\n")
  
  flog.info("Loaded portfolio with %s assets",ncol(p))
  out <- p[idx.inf:idx.sup, ]
  class(out) <- c(class(out), 'returns')

  if (is.null(rownames(out))) rownames(out) <- format(index(out), "%Y-%m-%d")
  out
}


