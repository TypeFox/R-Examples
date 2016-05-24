print.summary.heckitrob <-
function(x, ...)
{
  cat("------------------------------------------------------------- \n")
  cat("Robust 2-step Heckman / heckit M-estimation \n")
  cat("Probit selection equation: \n")
  names(x$coefficients$selection)=c("Estimate","Std.Error","t-value","p-value"," ")
  print(x$coefficients$selection)
  cat("Outcome equation: \n")
  names(x$coefficients$outcome)=c("Estimate","Std.Error","t-value","p-value"," ")
  print(x$coefficients$outcome)
  cat("---\n")
  cat("Signif. codes: 0 `***' 0.001 `**' 0.01 `*' 0.05 `.' 0.1 ` ' 1 \n")
  cat("-------------------------------------------------------------")
}
