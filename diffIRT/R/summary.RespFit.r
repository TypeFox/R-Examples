summary.RespFit=function(object,...){
  cat("\n")
  cat(if(toupper(object$model)==2) "Q-diffusion" else "D-diffusion","Model Fit of Responses\n")
  cat("---------------------------------\n")
  cat("Maydeu-Olivares & Joe Test of Order 2\n\n")

  cat("Univariate Statistics\n")
  print(object$Z)
  cat("\n\n")
 cat("Overall Test Statistic\n")
  cat("Mr =",round(object$Mr,3),"  df=",object$df,"  p = ",round(pchisq(object$Mr,object$df,lower.tail=F),3),"\n")
}

