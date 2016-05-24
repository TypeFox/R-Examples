mixpanelGetRetention <- function(
  account,
  event,                  # If missing, all events are included.
  from,
  to=from,
  unit="day",             # Defaults to 'day'.
  retentionType="birth",  # 'birth' (=First time) or 'compounded' (=Recurring). Defaults to 'birth'.
  intervalCount=15,       # Number of intervals per cohort to return.
  percentages=TRUE,       # Output as counts or percentages?
  ...                     # Additional arguments to Mixpanel API beside of <event>, <date_from>, <date_to>
                          # >> born_event='AppInstall'   # Needed for retention type 'birth'!!!
                          # >> born_where='properties["VersionString"]=="1.1"'
) {
  args = list(...)
  if (!missing(event))
    args$event = event
  
  args$from_date = createDateSequence(from)
  args$to_date = createDateSequence(to)
  args$unit = unit
  args$retention_type = retentionType
  args$interval_count = intervalCount
  
  data = mixpanelGetData(account, "retention/", args, data=TRUE)
  
  ## Returns a list of (count, retention[])-pairs. Retention arrays are of different lengths. 
  data = jsonlite::fromJSON(data)
  
  ## Put it into a matrix.
  counts = lapply(data, "[[", 2)
  retentions = lapply(data, "[[", 1)
  lens = unlist(lapply(retentions, length))
  maxLen = max(lens)
  dates = names(counts)
  
  data = matrix(NA, length(counts), 1 + maxLen, dimnames=list(dates, c("count", seq(maxLen) - 1)))
  for (i in seq(along=counts))
    data[i, 2:(1+lens[i])] = retentions[[i]]
  data[, 1] = unlist(counts)
  
  data = data[order(dates), , drop=FALSE]
  if (percentages)
    data[, -1] = data[, -1] / data[, 1] * 100
  data
}
