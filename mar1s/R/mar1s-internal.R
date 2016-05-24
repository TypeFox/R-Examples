`.split` <-
function(comp)
{
  if(NCOL(comp) > 1)
  {
    first <- comp[, 1]
    rest <- comp[, -1]
  }
  else
  {
    first <- comp
    rest <- NULL
  }

  return(list(first = first, rest = rest))
}

`.decomp` <-
function(object, start.time, xreg.absdata = NULL, init.absdata = NULL)
{
  if(!is.null(xreg.absdata))
  {
    xreg.cycl <- cycle(ts(xreg.absdata, start = start.time,
                          frequency = frequency(object$logseasonal)))
    xreg.logstoch <- log(xreg.absdata) -
      as.matrix(object$logseasonal)[xreg.cycl, -1]
  }
  else
    xreg.logstoch = NULL
  
  if(!is.null(init.absdata))
  {
    init.cycl <- cycle(ts(start = start.time -
                            deltat(object$logseasonal),
                          frequency = frequency(object$logseasonal)))
    init.logstoch <- log(init.absdata) -
      as.matrix(object$logseasonal)[init.cycl, ]
  }
  else
    init.logstoch = rep(0, NCOL(object$logseasonal))

  return(list(xreg.logstoch = xreg.logstoch,
              init.logstoch = init.logstoch))
}

`.mar1s.ts` <-
function(absdata, logdata, logstoch, logresid, names = NULL)
{
  colnames(absdata) <- names
  colnames(logdata) <- names
  colnames(logstoch)<- names
  result <- list(absdata = absdata, logdata = logdata,
                 logstoch = logstoch, logresid = logresid)
  class(result) <- 'mar1s.ts'
  return(result)
}
