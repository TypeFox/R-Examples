loop_data <-
function (data, id.x, id.y, nx) {
  # Get the x-th and y-th columns, only complete cases
  tmp <- data[, c(id.x, nx + id.y)]
  tmp <- tmp[complete.cases(tmp), ]
  
  x <- tmp [,1]
  y <- tmp [,2]
  names <- colnames(tmp)
  
  # Define the min/max, scope
  x.low  <- min(x)
  x.high <- max(x)
  y.low  <- min(y)
  y.high <- max(y)
  scope  <- (x.high - x.low) * (y.high - y.low)
  
  return (list(id.x=id.x, id.y=id.y, x=x, y=y, scope=scope, names=names, 
               x.low=x.low, x.high=x.high, y.low=y.low, y.high=y.high))
}