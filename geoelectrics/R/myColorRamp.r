#' Maps color to resistivity value
#' 
#' Maps color to (resistivity) values. 
#' A minimum and maximum value can be specified. 
#'   
#' @param col Character vector of colors.
#' @param values Numeric vector of values.
#' @param minData Minimum value (default min(values)). 
#' All smaller values will assigned to the first color in vector col.
#' @param maxData Maximum value (default max(values)).
#' All higher values will assigned to the last color in vector col.
#' @export
myColorRamp <- function(col, values, minData=min(values), maxData=max(values)) {
  v <- (values - minData) / diff(range(minData, maxData)) # same colors for all profiles
  v[v<0] <- 0
  v[v>1] <- 1
  x <- colorRamp(col)(v)
  rgb(x[,1], x[,2], x[,3], maxColorValue = 255) 
} 