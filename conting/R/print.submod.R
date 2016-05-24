print.submod <-
function(x,...,digits = max(3, getOption("digits") - 3)){

blob<-round(x$post_prob,digits=digits)

cat("Posterior model probability = ",blob,"\n")
cat("\n")

df<-data.frame(post_mean=x$post_mean,
post_var=x$post_var,lower_lim=x$lower,upper_lim=x$upper)
nam<-x$term
row.names(df)<-nam

lev<-100*x$prob.level

cat("Posterior summary statistics of log-linear parameters:\n")
print(df,digits=digits)
cat("NB: lower_lim and upper_lim refer to the lower and upper values of the\n")
cat(lev,"% highest posterior density intervals, respectively\n")
cat("\n")
if(!is.null(x$meanTOT)){
cat("Posterior mean of total population size =",round(x$meanTOT,digits),"\n")
cat(lev,"% highest posterior density interval for total population size = (",round(x$int,digits),") \n")
cat("\n")}

statistic2<-c("X2","deviance","Freeman-Tukey")[x$statnum==(1:3)]

cat("Under the",statistic2,"statistic \n")
cat("\n")
cat("Summary statistics for T_pred \n")
print(round(summary(x$Tpred),digits=digits))
cat("\n")
cat("Summary statistics for T_obs \n")
print(round(summary(x$Tobs),digits=digits))
cat("\n")
cat("Bayesian p-value = ",round(x$pval,digits=digits),"\n")

}
