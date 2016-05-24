rowModes <-
function(x) {
  # function for calculating the mode of all row values
  # analogous to rowMeans and rowMedians
  apply(x, MARGIN=2, Mode)
}
