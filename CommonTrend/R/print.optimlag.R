"print.optimlag" <-function(x,  ...){
  cat("Under the assumption of only one common trend," ,"\n",
      "the Lag order to maximize Johansen statistic(1988):","\n",sep="")
  print(x$Olag)
  cat("\n")
  cat("Corresponding Johansen Statistic:","\n")
  print(x$Olag.value)
  cat("\n")
  cat("A list for all the possible lag orders ", "\n")
  print(x$list.lags)
}