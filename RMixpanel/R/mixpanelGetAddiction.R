mixpanelGetAddiction <- function(
  account,
  event,                  # If missing, all events are included.
  from,
  to=from,
  unit="day",             # Defaults to 'day'.
  percentages=TRUE,       # Output as counts or percentages?
  addictionUnit="hour",   # Sub-unit for addiction.
  ...
) {
  args = list(...)
  if (!missing(event))
    args$event = event
  
  args$from_date = createDateSequence(from)
  args$to_date = createDateSequence(to)
  args$unit = unit
  args$addiction_unit = addictionUnit
  
  data = mixpanelGetData(account, "retention/addiction", args, data=TRUE)
  data = jsonlite::fromJSON(data)$data
  
  dates = names(data)
  data = matrix(unlist(data), length(data), byrow=TRUE, dimnames=list(dates, seq(along=data[[1]])))
  data = data[order(dates), , drop=FALSE]
  colnames(data)[1] = "count"
  
  if (percentages)
    data[, -1] = data[, -1] / data[, 1] * 100
  data
}
