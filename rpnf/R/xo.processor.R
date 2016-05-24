#' Determine the XO development of a given time series.
#' 
#' This is the main PNF-Workhorse, which transforms a given time series into X and O's.
#' Furthermore it provides already simple Buy/Sell signals, based on checking if the second last colum.
#' @param high
#' @param low
#' @param date
#' @param reversal
#' @param boxsize A single numeric value, indicating the boxsize to be considered.
#' @param log TRUE, if logarithmic scales should be used.
#' @keywords internal
xo.processor <- function(high,low=high, date, reversal=3L, boxsize=1, log=FALSE) {
  if (!(is.numeric(high) & length(high)>=1)) {
    stop("Argument high has to be numeric and at least of length 1!")
  }
  if (!(is.numeric(low) & length(low)>=1)) {
    stop("Argument low has to be numeric and at least of length 1!")
  }
  # FIXME this is not working yet
  # if (!(is.numeric.Date(date) & length(date)>=1)) {
  #   stop("Argument date has to be a Date object and has at least of length 1!")
  # }
  if (!(length(high)==length(low))) {
    stop("Arguments high and low have to have the same length!")
  }
  if (!(length(high)==length(date))) {
    stop("Arguments high and date have to have the same length!")
  }
  if (!(is.numeric(reversal) & length(reversal)==1)) {
    stop("Argument reversal has to be integer, and of length 1 exactly!")
  }
  if (!(is.integer(reversal) & min(reversal)>0)) {
    stop("Argument reversal has to be integer greater than zero!")
  }
  if (!(is.numeric(boxsize) & min(boxsize)>0)) {
    stop("Argument boxsize has to be numeric greater than zero, and of length 1 exactly!")
  }
  if (!(is.logical(log) & length(log)==1)) {
    stop("Argument log has to be numeric, and of length 1 exactly!")
  }
  
  # add additional columns
  status.xo <- rep(NA,length.out=length(high))  # current state of X/O signal
  boxnumber <- rep(NA,length.out=length(high))  # current integer box number (used for internal processes)
  status.bs <- rep(NA,length.out=length(high))  # current state of buy/sell signal
  nextX <- rep(NA,length.out=length(high))      # current nextX
  lastNextX <- rep(NA,length.out=length(high))  # remember last known nextX, necessary to identify buy/sell changes
  nextO <- rep(NA,length.out=length(high))      # current nextO
  lastNextO <- rep(NA,length.out=length(high))  # remember last known nextO, necessary to identify buy/sell changes
  column <- rep(NA,length.out=length(high))     # counter for current P&F column
  
  # TODO improve initialization
  status.xo[1] = "X"
  boxnumber[1] = quote2box(quote=high[1],boxsize=boxsize,log=log)
  status.bs[1] = "Buy"
  nextX[1] = high[1]
  lastNextX[1] = high[1]
  nextO[1] = low[1]
  lastNextO[1] = low[1]
  column[1] = 1
  
  # pnfprocessor
  for (i in 2:length(high)) {
    if (status.xo[i-1] == "X") {
      # we are in X-mode
      if (high[i] > nextX[i-1]) {
        # we made a new X
        status.xo[i] <- "X"
        boxnumber[i] = quote2box(quote=high[i],boxsize=boxsize,log=log)
        nextX[i] <- nextBox(high[i],"X", boxsize=boxsize,log=log)
        lastNextX[i] <- lastNextX[i-1]
        nextO[i] <- nextReversal(quote=high[i],status="X",reversal=reversal,boxsize=boxsize,log=log)
        lastNextO[i] <- lastNextO[i-1]
        column[i] = column[i-1]
        if (high[i] > lastNextX[i-1])
          status.bs[i] <- "Buy"
        else
          status.bs[i] <- status.bs[i-1]
      } else if (low[i] < nextO[i-1]) {
        # we made a reversal to O
        status.xo[i] <- "O"
        boxnumber[i] = quote2box(quote=low[i],boxsize=boxsize,log=log)
        nextX[i] <- nextReversal(low[i],"O",reversal=reversal,boxsize=boxsize,log=log)
        lastNextX[i] <- nextX[i-1]
        nextO[i] <- nextBox(low[i],"O", boxsize=boxsize,log=log)
        lastNextO[i] <- lastNextO[i-1]
        column[i] = column[i-1]+1
        if (low[i] < lastNextO[i-1])
          status.bs[i] <- "Sell"
        else
          status.bs[i] <- status.bs[i-1]
      } else {
        # nothing new happened
        status.xo[i] <- status.xo[i-1]
        boxnumber[i] = quote2box(quote=high[i],boxsize=boxsize,log=log)
        status.bs[i] <- status.bs[i-1]
        nextX[i] <- nextX[i-1]
        lastNextX[i] <- lastNextX[i-1]
        nextO[i] <- nextO[i-1]
        lastNextO[i] <- lastNextO[i-1]
        column[i] = column[i-1]
      }
    } else {
      # we are in O-mode
      if (low[i] < nextO[i-1]) {
        # we made a new O
        status.xo[i] <- "O"
        boxnumber[i] = quote2box(quote=low[i],boxsize=boxsize,log=log)
        nextO[i] <- nextBox(low[i],"O", boxsize=boxsize,log=log)
        lastNextO[i] <- lastNextO[i-1]
        nextX[i] <- nextReversal(low[i],"O",reversal=reversal,boxsize=boxsize,log=log)
        lastNextX[i] <- lastNextX[i-1]
        column[i] = column[i-1]
        if (low[i] < lastNextO[i-1])
          status.bs[i] <- "Sell"
        else
          status.bs[i] <- status.bs[i-1]        
      } else if (high[i] > nextX[i-1]) {
        # we made a reversal to X
        status.xo[i] <- "X"
        boxnumber[i] = quote2box(quote=high[i],boxsize=boxsize,log=log)
        nextO[i] <- nextReversal(high[i],"X",reversal=reversal,boxsize=boxsize,log=log)
        lastNextO[i] <- nextO[i-1]
        nextX[i] <- nextBox(high[i],"X", boxsize=boxsize,log=log)
        lastNextX[i] <- lastNextX[i-1]
        column[i] = column[i-1]+1
        if (high[i] > lastNextX[i-1])
          status.bs[i] <- "Buy"
        else
          status.bs[i] <- status.bs[i-1]
      } else {
        # nothing new happened
        status.xo[i] <- status.xo[i-1]
        boxnumber[i] = quote2box(quote=low[i],boxsize=boxsize,log=log)
        status.bs[i] <- status.bs[i-1]
        nextX[i] <- nextX[i-1]
        lastNextX[i] <- lastNextX[i-1]
        nextO[i] <- nextO[i-1]
        lastNextO[i] <- lastNextO[i-1]
        column[i] = column[i-1]
      }
    }
  }
  return (data.frame(date,high,low,boxnumber,column,status.xo,nextX,nextO,status.bs,lastNextX,lastNextO))
}
