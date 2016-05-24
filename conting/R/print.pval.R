print.pval <-
function(x, digits = max(3, getOption("digits") - 3),...){

statistic2<-c("X2","deviance","Freeman-Tukey")[x$statnum==(1:3)]

cat("Under the",statistic2,"statistic \n")
cat("\n")
cat("Summary statistics for T_pred \n")
print(round(summary(x$Tpred),digits=digits))
cat("\n")
cat("Summary statistics for T_obs \n")
print(round(summary(x$Tobs),digits=digits))
cat("\n")
cat("Bayesian p-value = ",round(x$pval,digits=digits),"\n")}
