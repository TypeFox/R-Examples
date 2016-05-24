### Gumbel distributions: borrowed from VGAM.
### amended to make dgumb(-Inf) = 0

dgumbel <- function(x, location = 0, scale = 1, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  zedd <- (x - location) / scale
  logdensity <- -zedd - exp(-zedd) - log(scale)
  logdensity[is.infinite(x)] <- -Inf
                           # works for both +Inf and -Inf
  if(log.arg)
     logdensity
  else {
     exp(logdensity)
  }
}

pgumbel <- function(q, location = 0, scale = 1) {
  answer <- exp(-exp(-(q - location) / scale))
  answer[scale <= 0] <- NaN
  answer
}
