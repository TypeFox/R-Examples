p_above <- 
function (loop.data, slope, intercept) {
  if (slope < 0) {
    return( NaN )
  }

  x <- loop.data$x
  y <- loop.data$y
  
  y.c <- slope * x + intercept
  
  # vertical difference between points on ceiling and observation for given x
  y.diff <- y - y.c
  
  for (i in 1:length(x)) {
    y.diff[i] <- max(0, y.diff[i])
  }
  
  # TODO Magic number
  return( sum(y.diff > 0.0000001) )
}