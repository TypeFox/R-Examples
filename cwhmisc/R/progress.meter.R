progress.meter <- function(i) {
  if      (i==0) cat("\n      0")
  else if (i %% 50 == 0) cat("\n",formatFix(i,0,6))
  else if (i %% 10 == 0) cat((i %/% 10) %% 10)
  else if (i %%  5 == 0) cat("+")
  else cat(".");
#  else cat(".");
  invisible(NULL)
}
