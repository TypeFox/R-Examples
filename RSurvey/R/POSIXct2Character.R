# Convert objects from POSIXct to character class

POSIXct2Character <- function(x, fmt="%Y-%m-%d %H:%M:%OS3") {
  pos <- gregexpr("%OS[[:digit:]]+", fmt)[[1]]
  if (pos > 0) {
    pos <- pos + c(3L,  attr(pos, "match.length"))
    dec.digits <- as.integer(substr(fmt, pos[1], pos[2]))
    x <- as.POSIXlt(x, tz=attr(x, "tzone"))
    x$sec <- round(x$sec, dec.digits) + 10^(-dec.digits - 1L)
  }
  return(format(x, format=fmt))
}
