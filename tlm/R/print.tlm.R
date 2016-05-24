print.tlm <-
function(x, ...)
 {
  printPreamble(x)
  print(x$model, ...)
  cat("\n")
 }
