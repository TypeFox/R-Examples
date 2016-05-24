readSudoku <- function(fn, map=c(1:9,letters)) {
  z <- scan(fn, "", quiet=TRUE)
  z <- do.call(rbind, strsplit(z, ""))
  matrix(match(z, map, 0), nrow(z))
}
