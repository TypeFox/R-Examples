Mode <-
function(x) {
  # function for calculating the mode
  ux = unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
