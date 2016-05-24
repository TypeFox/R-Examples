writeSudoku <- function(z, fn) {
  z[z==0] <- "-"
  cat(t(z), file=fn, sep=c(rep.int("", ncol(z)-1), "\n"))        # From "write"
}
