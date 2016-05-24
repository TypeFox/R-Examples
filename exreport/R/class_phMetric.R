is.phMetric <- function(x) {
  is(x, "phMetric")
}

print.phMetric <- function (x, ...) {
  print(x$data)
}

#Anonymous constructor
.phMetric <- function(data, format, decreasingOrder)
{
  metric <- list(
    "data"            = data,
    "format"          = format,
    "decreasingOrder" = decreasingOrder
  )
  
  class(metric) <- c("phMetric")
  
  metric
}

