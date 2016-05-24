fxtrim <- function(x, lo = 400, up = 10000)
{
  # convert latencies < 400 to NAs
  x[x < lo] <- NA
  # convert latencies > 10000 to NAs
  x[x > up] <- NA
  x
}
