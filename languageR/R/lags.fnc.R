lags.fnc = function(dat, time="Trial", group="Subject", depvar="RT", lag=1) {

  dat = dat[order(dat[,group],dat[, time]),]
  means = tapply(dat[,depvar], dat[,group], mean)

  f.fnc = function(dfr) {
    current = dfr[,depvar]
    previous = c(rep(as.numeric(means[dfr[1,group]]),lag), 
                 current[1:(length(current)-lag)])
    return(previous)
  }
  prev = by(dat, dat[,group], f.fnc, simplify=TRUE)
  return(unlist(prev))

}

