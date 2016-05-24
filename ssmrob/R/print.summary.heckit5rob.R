print.summary.heckit5rob <-
function(x, ...)
{
  cat("------------------------------------------------------------- \n")
  cat("Robust 2-step Heckman / heckit M-estimation \n")
  cat("Probit selection equation: \n")
  names(x$coefficients$selection)=c("Estimate","Std.Error","t-value","p-value"," ")
  print(x$coefficients$selection)
  cat("Outcome equation, regime 1: \n")
  names(x$coefficients$outcome1)=c("Estimate","Std.Error","t-value","p-value"," ")
  print(x$coefficients$outcome1)
  cat("Outcome equation, regime 2: \n")
  names(x$coefficients$outcome2)=c("Estimate","Std.Error","t-value","p-value"," ")
  print(x$coefficients$outcome2)
  cat("---\n")
  cat("Signif. codes: 0 `***' 0.001 `**' 0.01 `*' 0.05 `.' 0.1 ` ' 1 \n")
  cat("-------------------------------------------------------------")
}
