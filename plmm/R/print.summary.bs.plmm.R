print.summary.bs.plmm <-
function(x, ...){
  len<-length(x[[1]])
  coef_mat<-cbind(x$mean[-len], x$sd[-len])
  colnames(coef_mat)<-c("mean", "sd")
  cat("Coefficients:\n")
  print(coef_mat); cat("\n")
  cat("quantiles\n")
  print(t(x$quantiles[,-len])); cat("\n\n")
  
  VC_mat<-cbind(x$mean[len], x$sd[len])
  colnames(VC_mat)<-c("mean", "sd")
  cat("Variance component:\n")
  print(VC_mat); cat("\n")
  cat("quantiles\n")
  print(x$quantiles[,len]); cat("\n")
}
