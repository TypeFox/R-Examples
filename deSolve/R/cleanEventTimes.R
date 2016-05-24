## find nearest event for each time step
nearestEvent <- function(times, eventtimes) {
  eventtimes  <- unique(eventtimes) # remove double events first
  ## sorting does not cost much if already sorted
  times <- sort(times)
  eventtimes <- sort(eventtimes)
  ## find index of events where time is between
  inearest <- findInterval(times, eventtimes)
  ## special care for smallest and biggest element
  lower <- eventtimes[pmax(inearest, 1)]
  upper <- eventtimes[pmin(inearest + 1, length(eventtimes))]
  nearest <- ifelse(times - lower < upper - times, lower, upper)
  return(nearest)
}

## remove times that are numerically "too close" to an event
cleanEventTimes <- function(times, eventtimes, eps = .Machine$double.eps  * 10) {
  ## sorting does not cost much if already sorted
  ## sort times to ensure match of returned "nearest" value
  times <- sort(times)
  nearest <- nearestEvent(times, eventtimes)
  ## use bigger of the two numbers
  div <- pmax(times, nearest)
  ## special handling of zero
  div <- ifelse(div == 0, 1, div)
  reldiff <- abs(times - nearest) / div
  tooClose <- reldiff < eps
  times[!tooClose]
}
