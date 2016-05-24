`fit.mar1s` <-
function(x, xreg = NULL, seasonal.fun = seasonal.smooth, ...)
{
  absdata <- cbind(x, xreg)
  names <- colnames(absdata)
  
  logdata <- log(absdata)
  
  logseasonal <-
    do.call(cbind, lapply(logdata, seasonal.fun, ...))
  logstoch <- logdata -
    as.matrix(logseasonal)[cycle(logdata), ]
  y <- .split(logstoch)

  logstoch.ar1 <- arima(nan2na(y$first), xreg = nan2na(y$rest),
                        order = c(1, 0, 0), include.mean = FALSE)
  
  logresid.sd <- sd(resid(logstoch.ar1), TRUE)
  
  result <- list(logseasonal = logseasonal,
                 logstoch.ar1= logstoch.ar1,
                 logresid.sd = logresid.sd,
                 decomposed = .mar1s.ts(absdata, logdata, logstoch,
                   resid(logstoch.ar1), names))
  class(result) <- 'mar1s'
  return(result)
}
