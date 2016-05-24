print.gresponse.dmm <-
function(x, ...)
# print.gresponse.dmm() - print a gresponse.dmm object  (brief output)
{
  cat("Call:\n")
  print(x$call)
  cat("\nPredicted response to selection using component(s) ",x$effects,"\n")

  cat("\nGenetic selection differentials achieved by given psd:\n\n")
  cat("Overall:\n")
  print(x$overall,digits=x$digits)
  cat("\n")

}
