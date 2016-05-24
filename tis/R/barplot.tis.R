barplot.tis <- function(height, ...){
  ## If the first argument given to barplot() is a tis, this method is called
  ## There may be more than one series argument given, so the first thing to
  ## do is pull out the tis args.
  args <- c(list(height = height), list(...))
  tisArgs <- sapply(args, is.tis)
  series <- args[tisArgs]
  nonSeriesArgs <- args[!tisArgs]

  ncols <- sapply(series, NCOL)
  if(any(ncols != ncols[1])) stop("series have different column dimensions")
  zHeight <- t(do.call(rbind, lapply(lapply(series, stripTis), as.matrix)))
  dTimes <- function(x) diff(time(jul(c(start(x) - 1, ti(x)))))
  zWidth <- as.vector(unlist(sapply(series, dTimes)))

  beside <- is.logical(nonSeriesArgs$beside) && nonSeriesArgs$beside

  if(beside){
    ## if beside was set, just make each group of bars plus its space have width 1
    zWidth <- mean(zWidth)/(ncols[1] + 0.5)
    space <- c(0, 0.5)
  } else {
    ## set space argument
    space <- nonSeriesArgs$space
    if(is.null(space)) space <- sqrt(zWidth) * 0.2
    ## iterate a few times to set space and zWidth such that the width of a
    ## bar plus it's space equals the difference in times so that the
    ## bars are spaced correctly on the time axis.
    space <- rep(space, length.out = length(zWidth))
    zw <- zWidth
    for(i in 1:10){
      avg <- mean(zw, na.rm = T)
      zw <- zWidth - space * avg
    }
    zWidth <- zw
  }
  nonSeriesArgs$space <- space

  x.offset <- min(time(asTi(sapply(series, start)) - 1)) + zWidth[1]
  add <- is.logical(nonSeriesArgs$add) && nonSeriesArgs$add
  if(add){
    zArgs <- updateList(list(height = zHeight,
                             width  = zWidth,
                             x.offset  = x.offset,
                             axes = F),
                        nonSeriesArgs)
  } else {
    zArgs <- updateList(list(height = zHeight,
                             width  = zWidth,
                             x.offset  = x.offset,
                             xaxs = "r"),
                        nonSeriesArgs)
  }
  boink <- do.call(barplot2, zArgs)
  if(!add){
    max.x <- max(time(asTi(sapply(series, end))))
    Axis(x = c(x.offset, max.x), side = 1)
  }
  return(boink)
}
