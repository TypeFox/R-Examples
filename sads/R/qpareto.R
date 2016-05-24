qpareto <- function (p, shape, scale = min(p), lower.tail = TRUE, log.p = FALSE) {
  shape[shape <= 0] <- NaN
  scale[scale <= 0] <- NaN
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  scale/(1 - p)^(1/shape)
}
