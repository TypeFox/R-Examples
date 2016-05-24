rowMedians <-
function(x) {
  # function for calculating the median of all row values
  # analogous to rowSums and rowMeans
  apply(x, MARGIN=2, median, na.rm=T)
}
