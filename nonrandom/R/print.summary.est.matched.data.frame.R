
print.summary.est.matched.data.frame <- function(x,
                                                 ...)
{
  cat("\n         Summary for effect estimation \n\n")
  cat("\n Treatment/exposure:", x$treat)
  cat("\n Outcome:", x$outcome)
  cat("\n Effect measure:", x$meas, "\n")

  cat("\n Effect estimates:\n\n")
  print(x$eff.tab)
  cat("\n\n")
  
}
