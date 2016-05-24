`decompose.mar1s` <-
function(object, absdata, start.time = head(time(absdata), 1),
         init.absdata = rep(NA, NCOL(absdata)))
{
  absdata <- ts(absdata, start = start.time,
                frequency = frequency(object$logseasonal))

  logdata <- log(absdata)
  logstoch <- logdata -
    as.matrix(object$logseasonal)[cycle(logdata), ]
  
  arcoef <- head(coef(object$logstoch.ar1), 1)
  xregcoef <- tail(coef(object$logstoch.ar1), -1)  
  init.logstoch <- .decomp(object, start.time, NULL,
                           rep(NA, NCOL(absdata)))$init.logstoch
  y <- .split(logstoch)

  logresid <- decompose.ar1(arcoef, y$first, head(init.logstoch, 1),
                            xregcoef, y$rest, tail(init.logstoch, -1))

  return(.mar1s.ts(absdata, logdata, logstoch, logresid,
                   colnames(object$logseasonal)))
}
