printSudoku <- function(z) {
  z[z==0] <- " "
  for (r in 0:9) {
    if (r > 0) cat("  |", z[r,1:3], "|", z[r,4:6], "|", z[r,7:9], "|\n")
    if (r %% 3 == 0) cat("  +-------+-------+-------+\n")
  }
}
