fxwins <- function(x, lo = 300, up = 3000)
{
  # convert latencies < 300 to 300 and >3000 to 3000
  x[x < lo] <- lo
  x[x > up] <- up
  x
}