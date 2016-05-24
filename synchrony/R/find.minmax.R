find.minmax <- function (timeseries) {
  
  if (NCOL(timeseries)==1) {
    timeseries=cbind(1:NROW(timeseries), timeseries)
  }
  ## Find local maxima
  max.index <- find.minmax.aux (timeseries)
  ## Find local minima
  min.index <- find.minmax.aux (timeseries, mins=TRUE)
  
  mins=data.frame(index=min.index, time=timeseries[min.index,1], val=timeseries[min.index,2])
  maxs=data.frame(index=max.index, time=timeseries[max.index,1], val=timeseries[max.index,2])
  
  if (NCOL(mins)==1)
    mins=t(mins)
  if (NCOL(maxs)==1)
    maxs=t(maxs)
  
  rownames(mins)=1:NROW(mins)
  rownames(maxs)=1:NROW(maxs)
  results=list(mins=mins, maxs=maxs)
  class(results)="minmax"
  return (results)
}

find.minmax.aux <- function (timeseries, mins=FALSE) {
  if (mins)
    ind <- diff(c(Inf, timeseries[,2])) < 0
  else
    ind <- diff(c(-Inf, timeseries[,2])) > 0
  ind <- cumsum(rle(ind)$lengths)
  ind <- ind[seq.int(1, length(ind), 2)]
  ## First and last point cannot be local min/max
  ind = ind[!(ind %in% c(1, NROW(timeseries)))]
  return (ind)
}
