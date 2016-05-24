
print.summary.est.stratified.pscore <- function(x,
                                                ...)
{
  cat("\n         Summary for effect estimation \n\n")
  cat("\n Treatment/exposure:", x$treat)
  cat("\n Outcome:", x$outcome)
  cat("\n Effect measure:", x$meas, "\n")

  cat("\n Effect estimates:\n\n")
  print(x$eff.tab)

  cat("\n Stratum-specific parameter estimates:\n\n")
  print(format(x$str))

  cat("\n Stratum-specific adjusted parameter estimates:\n")
  cat(x$adjust)

  cat("\n Stratum-specific weights:\n")
  cat(x$weights)
  cat("\n\n")
  
}
