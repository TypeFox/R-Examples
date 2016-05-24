print.flimboot <-
function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nBootstrap flim-objects  are in $samples\n")
  cat("\nSee also ?plot.flimboot and ?flimSD")
}
