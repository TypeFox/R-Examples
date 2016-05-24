Ti_lander <- function(a) {
  # Given thrust a, find the optimal time Ti to fire the rocket
  # We assume the optimal time is less than Tmax
  Tmax <- 100
  g <- function(Ti) {
    l.out <- lander(Ti, a)
    return(100*min(l.out[1], 0) - max(l.out[2], 0))
  }
  Ti <- optimize(g, c(0, Tmax), maximum=T)$maximum
  return(Ti)
}
