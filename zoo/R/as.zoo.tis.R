as.zoo.tis <- function(x, class = "ti", ...) {
  if (class == "ti") {
    as.zoo(as.zooreg(x, class = "ti", ...))
  } else if (class == "numeric") {
    zoo(tis::stripTis(x), time(tis::ti(x), offset = 0))
  } else {
    asFun <- paste("as", class, sep = ".")
    zoo(tis::stripTis(x), do.call(asFun, list(tis::POSIXct(tis::ti(x), offset = 0, tz = "GMT"))), ...)
  }
}

as.zooreg.tis <- function(x, frequency = NULL, class = "ti", ...) {
  if (class == "ti") zooreg(tis::stripTis(x), start = start(x), ...)
    else as.zooreg(as.zoo(x, class = class, ...), frequency = frequency)
}

