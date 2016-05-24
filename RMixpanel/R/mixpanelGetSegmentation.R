mixpanelGetSegmentation <- function(
  account,
  event,                   # Must be included. E.g. '[Action] Video content view started'.
  from,
  to=from,
  unit="day",              # 
  type="unique",           # This can be "general", "unique", or "average".
  on='properties["$os"]',  # Array of up to 2 segmentation properties. An empty array returns un-segmented counts.
  action,                  # Could be "sum" or "average". If given, 1st property is aggregated by this function.
  verbose=TRUE,
  ...                      # Additional arguments to Mixpanel API.
) {
  args = list(...)
  args$event = event
  args$from_date = createDateSequence(from)
  args$to_date = createDateSequence(to)
  args$unit = unit
  args$type = type
  
  segmentDim = length(on)
  outDim = segmentDim + 1
  
  if (segmentDim > 2)
    stop("Up to 2 segmentation variables are handled by API.")
  
  if (!missing("action")) {
    args$action = action
    on[1] = paste('number(', on[1], ')', sep='') # Convert to numeric.
    outDim = outDim - 1  # Aggregation reduces dimension count.
  }
  
  if (segmentDim == 2) {
    args$inner = on[1]
    args$outer = on[2]
    data = mixpanelGetData(account, "segmentation/multiseg", args, data=TRUE, verbose=verbose)

  } else {
    args$on = on 
    data = mixpanelGetData(account, "segmentation/", args, data=TRUE, verbose=verbose)
  }
  
  ## API call.
  data = jsonlite::fromJSON(data)$data

  if (outDim == 3) {
    outerNames = names(data$values)
    innerNames = names(data$values[[1]])
    timeNames = names(data$values[[1]][[1]])
    
    kOuter = length(outerNames)
    kInner = length(innerNames)
    kTimes = length(timeNames)
    
    data = array(unlist(data$values), c(kTimes, kInner, kOuter), dimnames=list(timeNames, innerNames, outerNames))
    data[order(timeNames), , , drop=FALSE]
    
  } else { # outDim == 2 or 1.
    labels = names(data$values[[1]])
    n = length(labels)
    k = length(data[[2]])
    groups = names(data[[2]])
    
    data = matrix(unlist(data[[2]]), n, k, byrow=FALSE, dimnames=list(labels, groups))
    data[order(labels), , drop=FALSE]
  }
}
