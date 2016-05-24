"ts.update" <- 
function(x)
{
  ## convert an old-style ts, cts, its, or rts time series to a
  ## timeSeries or signalSeries object
  orig.time.zone <- timeDateOptions(time.zone = "GMT")
  on.exit(timeDateOptions(orig.time.zone))
  val <- if(is.ts(x)) {
    parms <- tsp(x)
    oldClass(x) <- NULL
    tsp(x) <- NULL
    units <- attributes(parms)
    if(!is.null(units))
      units <- units$units
    if(is.null(units))
      units <- character(0)
    x <- asSeriesData(x)
    names(parms) = c("start", "end", "deltat")
    if(!length(parms[["start"]]))
      parms[["start"]] <- 0
    signalSeries(data = x, units.position = units,
                 from = parms["start"], by = parms["deltat"])
  }
  else stop("Unknown type of time series")
  if(inherits(positions(val), "timeDate"))
    positions(val) <- timeZoneConvert(positions(val),
                                      orig.time.zone$time.zone)
  val
}

