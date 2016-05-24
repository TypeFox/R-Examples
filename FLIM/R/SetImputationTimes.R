SetImputationTimes <-
function(dataset, t.values) {
  if(is.null(t.values)) t.values <- sort(unique(dataset[, 2]))
  return(t.values)
}
