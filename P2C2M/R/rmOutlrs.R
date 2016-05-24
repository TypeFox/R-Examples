rmOutlrs <-
function(ind) {
  # function to remove outlier values
  outd = ind[!ind %in% boxplot.stats(ind)$out]
  return(outd)
}
