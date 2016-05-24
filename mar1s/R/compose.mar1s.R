`compose.mar1s` <-
function(object, loginnov, start.time = head(time(loginnov), 1),
         xreg.absdata = NULL, init.absdata = NULL)
{
  arcoef <- head(coef(object$logstoch.ar1), 1)
  xregcoef <- tail(coef(object$logstoch.ar1), -1)  
  d <- .decomp(object, start.time, xreg.absdata, init.absdata)

  y1 <- compose.ar1(arcoef, loginnov, head(d$init.logstoch, 1),
                    xregcoef, d$xreg.logstoch, tail(d$init.logstoch, -1))

  logstoch <- ts(cbind(y1, d$xreg.logstoch), start = start.time,
                 frequency = frequency(object$logseasonal))
  logdata <- logstoch +
    as.matrix(object$logseasonal)[cycle(logstoch), ]
  absdata <- exp(logdata)

  return(.mar1s.ts(absdata, logdata, logstoch, loginnov,
                   colnames(object$logseasonal)))
}
