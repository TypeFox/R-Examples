SlidingWindow <- function(FUN, data, window, step){
  total <- length(data)
  spots <- seq(from = 1, to = (total - window), by = step)
  result <- vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] <- match.fun(FUN)(data[spots[i]:(spots[i] + window - 1)])
  }
  return(result)
}