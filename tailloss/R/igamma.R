igamma <- function (s, u, log = FALSE) {
  robj <- pgamma(u, shape = s, rate = 1, log.p = TRUE) + lgamma(s)
  if (!log) robj <- exp(robj)
  robj
}
