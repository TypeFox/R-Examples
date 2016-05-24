
print.summary.relative.effect <- function(x,
                                          ...)
{
  cat("\n Outcome: ", x$resp, "\n")
  cat("\n Treatment: ", x$treat, "\n")
  cat("\n Covariates: ", x$sel, "\n")
  cat("\n Unadjusted treatment effect: ", x$unadj, "\n")
  cat("\n Results:\n")
  print(x$result)

}
