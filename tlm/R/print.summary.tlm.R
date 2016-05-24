print.summary.tlm <-
function(x, ...)
 {
  printPreamble(x)
  print(x$summary, ...)
  cat("\n")
 }
