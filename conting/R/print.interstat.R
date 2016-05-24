print.interstat <-
function(x,digits=max(3, getOption("digits") - 3),...){

df<-data.frame(post_prob=x$prob,post_mean=x$post_mean,
post_var=x$post_var,lower_lim=x$lower,upper_lim=x$upper)
nam<-x$term
row.names(df)<-nam

lev<-100*x$prob.level

cat("Posterior summary statistics of log-linear parameters:\n")
print(df,digits=digits)
cat("NB: lower_lim and upper_lim refer to the lower and upper values of the\n")
cat(lev,"% highest posterior density intervals, respectively\n")}
